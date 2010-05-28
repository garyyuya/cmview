package cmview.gmbp;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Iterator;
import java.util.Vector;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import javax.swing.BorderFactory;
import javax.swing.JPanel;

import owl.core.structure.features.SecStrucElement;
import owl.core.structure.graphs.RIGNbhood;
import owl.core.structure.graphs.RIGNode;

import owl.gmbp.CMPdb_nbhString_traces;
import owl.gmbp.CMPdb_sphoxel;
import owl.gmbp.CSVhandler;

import cmview.ContactMapPane;
import cmview.ScreenBuffer;
import cmview.Start;
import cmview.datasources.Model;
import edu.uci.ics.jung.graph.util.Pair;

public class ContactPane extends JPanel implements MouseListener, MouseMotionListener, ComponentListener{ //KeyListener
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	protected static final Dimension defaultDim = new Dimension(1700, 850);
	protected static final float g2dRatio = 0.5f; // H/W
	protected static final double defaultMinAllowedRat = -3;
	protected static final double defaultMaxAllowedRat = 1;
	protected final int cylindricalMapProj = 0;
	protected final int kavrayskiyMapProj = 1;
	protected final int azimuthalMapProj = 2;
	protected final int defaultProjType = kavrayskiyMapProj;
	protected static final int sphoxelHist = 0;
	protected static final int tracesHist = 0;

	protected static final Color helixSSTColor = Color.magenta; // color for helix residues
	protected static final Color sheetSSTColor = Color.yellow;	// color for helix residues
	protected static final Color otherSSTColor = Color.cyan;	// color for helix residues
	protected static final Color anySSTColor = Color.blue;		// color for helix residues
		
	/*--------------------------- member variables --------------------------*/		
	// underlying data
	private Model mod;
	private ContactMapPane cmPane;
	private ContactView contactView;
	
	private ContactStatusBar contStatBar;	
	
	// used for drawing
	private Dimension screenSize;			// current size of this component (contactView) on screen
	private Dimension g2dSize;         		// size of the effective Rectangle
											// available for drawing the sphoxel image
//	private double screenRatio;  			// ratio of screen size and contact
//											// map size = size of each contact
											// on screen
//	private int sphoxelSize;
	
	private Point mousePressedPos;   		// position where mouse where last
											// pressed, used for start of square
											// selection and for common
											// Neighbors
	private Point mouseDraggingPos;  		// current position of mouse
											// dragging, used for end point of
											// square selection
	private Point mousePos;             	// current position of mouse (being
											// updated by mouseMoved)
	private int lastMouseButtonPressed;		// mouse button which was pressed
											// last (being updated by
											// MousePressed)
//	private int currentRulerCoord;	 		// the residue number shown if
//											// showRulerSer=true
//	private int currentRulerMousePos;		// the current position of the mouse
//											// in the ruler
	private boolean dragging;     			// set to true while the user is
											// dragging (to display selection
											// rectangle)
	protected boolean mouseIn;				// true if the mouse is currently in
											// the contact map window (otherwise
											// supress crosshair)
//	private boolean showRulerCoord; 		// while true, current ruler
//											// coordinate are shown instead of
//											// usual coordinates
//	private boolean showRulerCrosshair;		// while true, ruler "crosshair" is
//											// being shown
	
	// buffers for triple buffering
	private ScreenBuffer screenBuffer;		// buffer containing the more or
											// less static background image
	
	// drawing colors (being set in the constructor)
	private Color backgroundColor;	  		// background color
	private Color squareSelColor;	  		// color of selection rectangle
	private Color crosshairColor;     		// color of crosshair	
	private Color selAngleRangeColor;       // color for selected rectangles
	private Color longitudeColor;			// color for longitudes
	private Color latitudeColor;			// color for latitudes
	
	// selections 
	private Vector<Pair<Double>> lampdaRanges;			// permanent list of currently selected lampda ranges
	private Vector<Pair<Double>> phiRanges;       // permanent list of currently selected phi ranges
	private Vector<Pair<Integer>> selContacts;         // permanent list of currently selected and referred contacts
	private Pair<Double> tmpLampdaRange;
	private Pair<Double> tmpPhiRange;
	private Pair<Integer> tmpSelContact;
	private Pair<Double> tmpAngleComb;
	private Pair<Double> origCoordinates;    // actual coordinates of origin (equator, zero meridian) in angle range x:[0:2*PI] y:[0:PI]
	private Pair<Double> centerOfProjection; // longitude and latitude of actual centre of projection (depends on origCoordinates)
	private Pair<Double> rightClickAngle;
	
	private RIGNode nodeI, nodeJ;
	
	// query data
	private char iRes='A', jRes='A';
	private char iSSType='H', jSSType='H';
	private boolean diffSStype=false;
//	private String iResType="Ala", jResType="Ala";
	private int iNum=0, jNum=0;
	private String nbhString, nbhStringL;
	private String origNBHString; 
	private String jAtom = "CA";
	private char[] nbhsRes;
	
	private String db = "bagler_all13p0_alledges";

	// Sphoxel-Data
	private CMPdb_sphoxel sphoxel;
	private double [][] ratios;
	private double [][][] bayesRatios;
	private double minRatio = 0;
	private double maxRatio = 0;
	private float minr=2.0f, maxr=12.8f;
	private float[] radiusThresholds = new float[] {2.0f, 5.6f, 9.2f, 12.8f};
	private boolean radiusRangesFixed = true;
	private String radiusPrefix = "rSR";
	
	private int [] histWholeAngleRange; // = new int[scale.length];
	private Vector<int[]> hist4SelAngleRange; // = new Vector<int[]>();
	private int [][] histTWholeAngleRange; // = new int[scale.length];
	private Vector<int[][]> histT4SelAngleRange; // = new Vector<int[]>();
	private Vector<double[]> minMaxAverT4SelAngleRange;
	private boolean showHist = true;
	private int histType = sphoxelHist;	
	private int chosenSelection;
	
	// NBHStraces-Data
	private CMPdb_nbhString_traces nbhsTraces;
	private Vector<float[]> nbhsNodes;
	private int[] numNodesPerLine;
	private int[] maxDistsJRes, minDistsJRes;
	private int numLines;
	
	// ---- variables for representation (drawing methods)
	private int numSteps = CMPdb_sphoxel.defaultNumSteps; //CMPdb_sphoxel.defaultNumSteps;  // change later via interface
	private float resol = CMPdb_sphoxel.defaultResol; //CMPdb_sphoxel.getDefaultResol();
	private final int border = 0; //15;
	private final int yBorderThres = 0; //45;
	private float pixelWidth = 5*36/this.numSteps; // change in dependence of numSteps
	private float pixelHeight = this.pixelWidth; //5*36/this.numSteps;	
	private float voxelsize = (float) (this.numSteps*this.pixelWidth/Math.PI); //pixelWidth;
	private final double voxelsizeFactor = 1.5;
	private double deltaRad = Math.PI/this.numSteps;
	private final double dL = 0.5;
	private double rSphere;
	private double maxDistPoints;

	private boolean removeOutliers = false; //true;
	private double minAllowedRat = defaultMinAllowedRat;
	private double maxAllowedRat = defaultMaxAllowedRat;
	private int chosenColourScale = ContactStatusBar.BLUERED;
	
	private boolean paintCentralResidue = true;
//	private boolean mapProjection = false;
	private int mapProjType = defaultProjType;
	private double deltaOffSetX, deltaOffSetXEnd, deltaOffSetXCenter; // = this.getScreenPosFromRad(0, 0).getFirst();
	
	private int xDim=0;
	private int yDim=0;
	
	public final char[] aas = new char[]{'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}; // possible Residues
//	private final String AAStr = new String(aas); 
	public final char[] sstypes = new char[]{'H','S','O','A'};
//	private final String SSTStr = new String(sstypes);
	
	/*----------------------------- constructors ----------------------------*/

	/**
	 * Create a new ContactPane.
	 * 
	 * @param mod
	 * @param cmPane
	 * @param view
	 */
	public ContactPane(Model mod, ContactMapPane cmPane, ContactView contactView){
		this.mod = mod;
		this.cmPane = cmPane;
		this.contactView = contactView;		
		addMouseListener(this);
		addMouseMotionListener(this);
		addComponentListener(this);
//		this.setFocusable(true);
//		addKeyListener(this);

		this.setOpaque(true); // make this component opaque
		this.setBorder(BorderFactory.createLineBorder(Color.black));
		
		// initialize data
		this.nbhsNodes = new Vector<float[]>();
		this.ratios = new double[0][0];
				
		try {
			calcParam();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		
		this.xDim = defaultDim.width;
//		this.yDim = defaultDim.height;
		// update pixel dimensions
		this.pixelWidth = (float)(this.xDim-2*this.border)/(float)(2*this.numSteps) ;
		this.pixelHeight =  this.pixelWidth;
		this.voxelsize = (float) ((float) (this.numSteps*this.pixelWidth)/Math.PI);
	
		this.g2dSize = defaultDim; // new Dimension(this.xDim,this.yDim);
		this.setSize(this.g2dSize);
		this.screenSize = this.contactView.getPreferredSize();
//		System.out.println("ContactPane screensize HxW: "+this.screenSize.height+"x"+this.screenSize.width);
//		this.screenSize = this.contactView.getScreenSize();
//		System.out.println("ContactPane getscreensize HxW: "+this.screenSize.height+"x"+this.screenSize.width);
//		this.contactView.setPreferredSize(this.screenSize);
//		this.contactView.setPreferredSize(new Dimension(this.screenSize.width+AngleRuler.STD_RULER_WIDTH, this.screenSize.height+AngleRuler.STD_RULER_WIDTH));

		this.rSphere = g2dSize.getHeight()/2;
		this.maxDistPoints = distPointsOnSphere(Math.PI, Math.PI/2, Math.PI, 0);
		
		this.mousePos = new Point();
		this.mousePressedPos = new Point();
		this.mouseDraggingPos = new Point();
		
		// set default colors
		this.backgroundColor = Color.white;
		this.squareSelColor = Color.gray;	
		this.crosshairColor = Color.green;
		this.selAngleRangeColor = Color.black;
		this.longitudeColor = Color.black;
		this.latitudeColor = this.longitudeColor;
		
		this.dragging = false;
		this.selContacts = new Vector<Pair<Integer>>();
		this.lampdaRanges = new Vector<Pair<Double>>();
		this.phiRanges = new Vector<Pair<Double>>();
		
		this.tmpLampdaRange = new Pair<Double>(0.0, 0.0);
		this.tmpPhiRange = new Pair<Double>(0.0, 0.0);
		this.tmpAngleComb = new Pair<Double>(0.0, 0.0);
		this.origCoordinates = new Pair<Double>(Math.PI, Math.PI/2);
//		this.centerOfProjection = new Pair<Double>(Math.PI-this.origCoordinates.getFirst(), Math.PI/2-this.origCoordinates.getSecond());
//		this.centerOfProjection = new Pair<Double>(this.origCoordinates.getFirst()-Math.PI, this.origCoordinates.getSecond()-Math.PI/2);
		this.centerOfProjection = this.origCoordinates;
		this.deltaOffSetX = this.getOffSet(0.0, 0.0); //this.getScreenPosFromRad(0, 0).getFirst();
		this.deltaOffSetXEnd = this.getOffSet(2*Math.PI, 0.0); //this.getScreenPosFromRad(2*Math.PI, 0.0).getFirst();
		this.deltaOffSetXCenter = this.getOffSet(Math.PI, 0.0);
		
//		setOutputSize(screenSize.width, screenSize.height);
//		System.out.println("outputsize= "+this.outputSize);
	}

	private void updateAngleRange(){
//		RIGEdge edge = this.mod.getGraph().getEdgeFromSerials(this.tmpSelContact.getFirst(), this.tmpSelContact.getSecond());		
//		edge.setPhiPsi(this.tmpLampdaRange.getFirst(), this.tmpLampdaRange.getSecond(), this.tmpPhiRange.getFirst(), this.tmpPhiRange.getSecond());		
		this.mod.getGmbp().setLampdaRanges(this.lampdaRanges);
		this.mod.getGmbp().setPhiRanges(this.phiRanges);
		this.mod.getGmbp().setSelContacts(this.selContacts);
	}
	
	public void updateQueryParam(int type){
		switch (type){
		case 0:
			calcSphoxelParam(); break;
		case 1:
			calcTracesParam(); break;
		}
	}
	
	private void calcParam() throws SQLException{
		calcSphoxelParam();
		this.sphoxel = new CMPdb_sphoxel(this.iRes, this.jRes, this.db);
//		setSphoxelParam(); // performed within calcSphoxel		
		calcSphoxel();

		calcTracesParam();	
		this.nbhsTraces = new CMPdb_nbhString_traces(this.nbhStringL, this.jAtom, this.db);
		setTracesParam();
		calcNbhsTraces();	
	}
	
	private void calcTracesParam(){
		RIGNbhood nbhood = this.mod.getGraph().getNbhood(nodeI);
		System.out.println("Edge type: "+this.mod.edgeType);
		this.jAtom = this.mod.edgeType.toUpperCase();
		this.nbhString = nbhood.getNbString();
		this.nbhStringL = "%";
		int count = 0;
//		int indexOfX = 0;
		this.nbhsRes = new char[this.nbhString.length()];
		for (int i=0; i<this.nbhString.length(); i++){
			this.nbhStringL += this.nbhString.charAt(i);
			this.nbhStringL += "%";
			this.nbhsRes[count] = this.nbhString.charAt(i);
//			if (this.nbhString.charAt(i) == 'x')
//				indexOfX = count;
			count++;
		}
		System.out.println(this.nbhString+"-->"+this.nbhStringL);	
		this.origNBHString = this.nbhString;
//		this.origNBHStringL = this.nbhStringL;
		
//		String s = "";
//		computeNBHStringCombinations(0, s);
	}
	
//	private void computeNBHStringCombinations(int index, String s){
//		String newS = s;
//		if (this.nbhsRes[index] != 'x'){
//			// this.nbhsRes[index] --> false
//			if (index+1 >= this.nbhsRes.length)
//				System.out.println("String: "+newS);
//			else
//				computeNBHStringCombinations(index+1, newS);			
//		}
//		// this.nbhsRes[index] --> true
//		if (index+1 >= this.nbhsRes.length)
//			System.out.println("String: "+newS);
//		else
//			computeNBHStringCombinations(index+1, newS+this.nbhsRes[index]);
//	}
//	
//	private String[] computeNBHStringCombinations(String[] comb, int length){
//		String[] newComb = new String[comb.length+1];
//		
//		return newComb;
//	}
	
	public void calcSphoxelParam(){
		// Get first shell neighbours of involved residues
		Pair<Integer> currentResPair = this.cmPane.getmousePos();
		calcSphoxelParam(currentResPair);
	}
		
	public void calcSphoxelParam(Pair<Integer> currentResPair){
		this.iNum = currentResPair.getFirst();
		this.jNum = currentResPair.getSecond();
//		System.out.println("first:"+this.iNum+"  second:"+this.jNum);
		
		// use pair to get iRes and jRes, isstype, nbhstring
		this.nodeI = this.mod.getNodeFromSerial(this.iNum); //this.mod.getGraph().getNodeFromSerial(this.iNum);
		nodeJ = this.mod.getNodeFromSerial(this.jNum);		
		this.iRes = this.nodeI.toString().charAt(this.nodeI.toString().length()-1);
		this.jRes = nodeJ.toString().charAt(nodeJ.toString().length()-1);
//		System.out.println("this.nodeI="+this.nodeI.toString()+"-->"+iRes+"  this.nodeI="+nodeJ.toString()+"-->"+jRes);
		
//		this.iResType = this.nodeI.getResidueType();
//		this.jResType = this.nodeJ.getResidueType();
//		System.out.println("iresType: "+this.iResType+"  jresType: "+this.jResType);

		// Definition of sstype and jatom
		SecStrucElement iSSelem = this.nodeI.getSecStrucElement();
//		type = this.nodeI.getSecStrucElement().getType(); 
		if (iSSelem == null){
			System.out.println("No SSelement!");
			this.diffSStype = false;
			this.iSSType = CMPdb_sphoxel.AnySStype;
		}
		else{
			if (iSSelem.isHelix())
				this.iSSType = SecStrucElement.HELIX;
			else if (iSSelem.isOther())
				this.iSSType = SecStrucElement.OTHER;
			else if (iSSelem.isStrand())
				this.iSSType = SecStrucElement.STRAND;
			else if (iSSelem.isTurn())
				this.iSSType = SecStrucElement.TURN;
			
			this.diffSStype = true;
		}
//		System.out.println("i secStrucElement: "+this.iSSType);
//		SecStrucElement jSSelem = nodeJ.getSecStrucElement();
//		if (jSSelem == null){
//			System.out.println("No JSSelement!");
//			this.jSSType = this.sphoxel.AnySStype;
//		}
//		else{
//			if (jSSelem.isHelix())
//				this.jSSType = SecStrucElement.HELIX;
//			else if (jSSelem.isOther())
//				this.jSSType = SecStrucElement.OTHER;
//			else if (jSSelem.isStrand())
//				this.jSSType = SecStrucElement.STRAND;
//			else if (jSSelem.isTurn())
//				this.jSSType = SecStrucElement.TURN;
//		}
		// Standard jSStype=any --> no differentiation made
		this.jSSType = CMPdb_sphoxel.AnySStype;
//		System.out.println("j secStrucElement: "+this.jSSType);	
	}
	
	private void setSphoxelParam(){
		this.sphoxel.setDiffSSType(this.diffSStype); // set to true if you want to differentiate ssType
		this.sphoxel.setISSType(this.iSSType);
		this.sphoxel.setJSSType(this.jSSType);
		this.sphoxel.setIRes(this.iRes);
		this.sphoxel.setJRes(this.jRes);
		this.sphoxel.setNumSteps(this.numSteps); // choose number of steps for resolution
		this.sphoxel.setMinr(this.minr);
		this.sphoxel.setMaxr(this.maxr);
//		this.ratios = new double [this.numSteps][2*this.numSteps];
	}
	
	private void calcSphoxel() throws SQLException{
		this.ratios = new double [this.numSteps][2*this.numSteps];
		
		if (this.radiusRangesFixed){

			if (this.minr==this.radiusThresholds[0] && this.maxr==this.radiusThresholds[1])
				this.radiusPrefix = CMPdb_sphoxel.radiusRanges[0];
			else if (this.minr==this.radiusThresholds[1] && this.maxr==this.radiusThresholds[2])
				this.radiusPrefix = CMPdb_sphoxel.radiusRanges[1];
			else if (this.minr==this.radiusThresholds[2] && this.maxr==this.radiusThresholds[3])
				this.radiusPrefix = CMPdb_sphoxel.radiusRanges[2];
			
			File dir1 = new File (".");		
			String fn = "/Users/vehlow/Documents/workspace/CMView"+Start.SPHOXEL_DIR; // 
//			fn = "/Users/vehlow/Documents/workspace/outputFiles/sphoxelBG/";
			try {
				fn = dir1.getCanonicalPath() + Start.SPHOXEL_DIR;
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			String archiveFN = fn + "SphoxelBGs.zip";
			ZipFile zipfile = null;
			try {
				zipfile = new ZipFile(archiveFN);
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			fn = "";
			
			fn = fn+"sphoxelBG_"+this.iRes+"-"+String.valueOf(this.iSSType).toLowerCase()+"_"+this.jRes+"-"
				+String.valueOf(CMPdb_sphoxel.AnySStype).toLowerCase()+"_"+this.radiusPrefix+".csv";
			System.out.println("Filename= "+fn);
//			fn = "/Users/vehlow/Documents/workspace/outputFiles/LogOddsScoresBayes_fromDB-bagler_all13p0_alledges_A-A_SStype-H_radius9.2-12.8_resol90.csv";
			
			ZipEntry zipentry = zipfile.getEntry(fn);
			
			CSVhandler csv = new CSVhandler();
			try {
//				this.bayesRatios = csv.readCSVfile3Ddouble(fn);
				this.bayesRatios = csv.readCSVfile3Ddouble(zipfile, zipentry);
				setNumSteps(this.bayesRatios.length);
				this.ratios = new double[this.bayesRatios.length][this.bayesRatios[0].length];
				for (int i=0; i<this.bayesRatios.length; i++){
					for (int j=0; j<this.bayesRatios[i].length; j++){
						this.ratios[i][j] = this.bayesRatios[i][j][0];
					}
				}
				for(int i=0;i<this.ratios.length;i++){ // dim for phi
					for(int j=0;j<this.ratios[i].length;j++){ // dim for lampda					                 
						if (this.ratios[i][j]<this.minRatio)
							this.minRatio = this.ratios[i][j];
						if (this.ratios[i][j]>this.maxRatio)
							this.maxRatio = this.ratios[i][j];
					}  
				}	
			} catch (NumberFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		else {
			setSphoxelParam();
							
			this.sphoxel.runBayesPreComp();
			this.ratios = this.sphoxel.getRatios();
			System.out.println("BayesRatios computed");		
//			this.bayesRatios = sphoxel.getBayesRatios();
			this.minRatio = this.sphoxel.getMinRatio();
			this.maxRatio = this.sphoxel.getMaxRatio();
		}
			
	}
	
	public void recalcSphoxel() throws SQLException{
		if (this.sphoxel.getNumSteps()!=this.numSteps || this.sphoxel.getMinr()!=this.minr || this.sphoxel.getMaxr()!=this.maxr 
				|| this.sphoxel.getIRes()!=this.iRes || this.sphoxel.getJRes()!=this.jRes 
				|| this.sphoxel.getISSType()!=this.iSSType || this.sphoxel.getJSSType()!=this.jSSType || this.sphoxel.getDiffSSType()!=this.diffSStype){
//			setSphoxelParam();
			calcSphoxel();
			
			updateScreenBuffer();
		}
	}
	
	private void setTracesParam(){
		nbhsTraces.setDiffSSType(this.diffSStype);
//		nbhsTraces.setDiffSSType(false);
		nbhsTraces.setSSType(this.iSSType);
		nbhsTraces.setNBHS(this.nbhStringL);
	}
	
	private void calcNbhsTraces() throws SQLException{
//		nbhsTraces.run();
		System.out.println("NbhsTraces extracted");
		nbhsNodes = nbhsTraces.getNBHSnodes();
		
		// FAKE
		CSVhandler csv = new CSVhandler();
		String sFileName = "/Users/vehlow/Documents/workspace/outputFiles/NBHSnodes_fromDB-bagler_cb8p0_alledges_nbhs-%C%P%x%H%G%_JAtom-CA.csv";
		if (sFileName!=null){
            System.out.println("Chosen path/file:" + sFileName);
		}
		else
            System.out.println("No path chosen!");
		try {
			this.nbhsNodes = csv.readCSVfileVector(sFileName);
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		this.nbhString = "CPxHG";
		this.nbhStringL = "%";
		int count = 0;
		this.nbhsRes = new char[this.nbhString.length()];
		for (int i=0; i<this.nbhString.length(); i++){
			this.nbhStringL += this.nbhString.charAt(i);
			this.nbhStringL += "%";
			this.nbhsRes[count] = this.nbhString.charAt(i);
			count++;
		}
		// END FAKE
		
		// compute nbhstringTraces		
		if (this.nbhsNodes.size()>0){
			System.out.println("this.nbhsNodes.size(): "+this.nbhsNodes.size());
			float[] node, nbNode;
			int numNodes = 1, minDist = 1, maxDist = 1;
			Vector<Integer> numNodesL = new Vector<Integer>();	
			Vector<Integer> minDists = new Vector<Integer>();	
			Vector<Integer> maxDists = new Vector<Integer>();	
			this.numLines = 1;
			node = (float[]) this.nbhsNodes.get(0);
			minDist = (int) node[2];
			for(int i=1; i<this.nbhsNodes.size(); i++){
				node = (float[]) this.nbhsNodes.get(i);	
				nbNode = (float[]) this.nbhsNodes.get(i-1);
//				System.out.println(node[0]+"-"+nbNode[0] +"  "+ node[1]+"-"+nbNode[1]);
				if (node[0]==nbNode[0] && node[1]==nbNode[1] && i>0){
					numNodes++;
				}
				else{
//				if (node[0]!=nbNode[0] || node[1]!=nbNode[1]){ //nodes not of same line
					minDists.add(minDist);// add first before update mindist
					minDist = (int) node[2];
					maxDist = (int) nbNode[2];// first update maxdist before adding
					maxDists.add(maxDist);
					numNodesL.add(numNodes);
					this.numLines++;
					numNodes = 1;
				}
			}
			numNodesL.add(numNodes);
			minDists.add(minDist);
			maxDists.add(maxDist);
			this.numNodesPerLine = new int[this.numLines];
			this.minDistsJRes = new int[this.numLines];
			this.maxDistsJRes = new int[this.numLines];
			System.out.println(this.numLines +" =? "+numNodesL.size());
			for(int i=0; i<numNodesL.size(); i++){
				this.numNodesPerLine[i] = (int) numNodesL.get(i);
				this.minDistsJRes[i] = (int) minDists.get(i);
				this.maxDistsJRes[i] = (int) maxDists.get(i);
//				System.out.print(this.numNodesPerLine[i]+"_"+this.minDistsJRes[i]+"_"+this.maxDistsJRes[i]+"\t");
			}
//			System.out.println();
		}	
		
	}
	
	public void recalcTraces(boolean perform) throws SQLException{
		if (perform){
			setTracesParam();
			calcNbhsTraces();
			
			updateScreenBuffer();			
		}
	}
	
	public void recalcTraces() throws SQLException{
		if (this.nbhsTraces.getDiffSSType()!=this.diffSStype || this.nbhsTraces.getSSType()!=this.iSSType){
			recalcTraces(true);
		}
	}
	
	public void calcHistogramms(){
		calcHistOfSphoxelMap();
		calcHistOfNBHSTraces();
	}
	
	public void calcHistOfNBHSTraces(){
//		ColorScaleView colView = this.contactView.getColorScaleView();
//		double [] scale = colView.getScaleValues();
		this.histTWholeAngleRange = new int[4][aas.length];
		this.histT4SelAngleRange = new Vector<int[][]>();
		this.minMaxAverT4SelAngleRange = new Vector<double[]>();
		int[][][] tempHist = new int[this.lampdaRanges.size()][4][aas.length];
		double[] averPhi = new double[this.phiRanges.size()];
		double[] averLampda = new double[this.lampdaRanges.size()];
		double[] minPhi = new double[this.lampdaRanges.size()];
		double[] maxPhi = new double[this.lampdaRanges.size()];
		double[] minLampda = new double[this.lampdaRanges.size()];
		double[] maxLampda = new double[this.lampdaRanges.size()];
		int[] foundNodes = new int[this.lampdaRanges.size()];
		
		float[] node;
		float phiRad, lampdaRad;
		double tMin, tMax, pMin, pMax;
		int jResID, jSSTypeID;
//		char jRes, jSSType;
		
		for(int i=0; i<this.nbhsNodes.size(); i++){			
			node = (float[]) this.nbhsNodes.get(i);
			phiRad = node[3];
			lampdaRad = (float) (node[4] + Math.PI);
			jResID = (int) node[5];
//			jRes = this.aas[jResID];
			jSSTypeID = (int) node[6];
//			jSSType = this.sstypes[jSSTypeID];
			
			this.histTWholeAngleRange[jSSTypeID][jResID] += 1; 
			this.histTWholeAngleRange[this.sstypes.length-1][jResID] += 1; 
			
			Iterator<Pair<Double>> itrP = this.lampdaRanges.iterator();
			Iterator<Pair<Double>> itrT = this.phiRanges.iterator();
			int index = 0;
			while (itrP.hasNext()){
				Pair<Double> lampda = (Pair<Double>) itrP.next();
				Pair<Double> phi = (Pair<Double>) itrT.next();
				tMin = phi.getFirst();
				tMax = phi.getSecond();
				pMin = lampda.getFirst();
				pMax = lampda.getSecond();
				
				if (phiRad>=tMin && phiRad<tMax && lampdaRad>=pMin && lampdaRad<pMax){
					tempHist[index][jSSTypeID][jResID] += 1; 
					tempHist[index][this.sstypes.length-1][jResID] += 1; 
					
					if (phiRad>maxPhi[index])
						maxPhi[index]=phiRad;
					if (lampdaRad>maxLampda[index])
						maxLampda[index]=lampdaRad;
					if (phiRad<minPhi[index] || minPhi[index]==0)
						minPhi[index]=phiRad;
					if (lampdaRad<minLampda[index] || minLampda[index]==0)
						minLampda[index]=lampdaRad;
					averPhi[index] += phiRad;
					averLampda[index] += lampdaRad;
					foundNodes[index] += 1;
				}
				
				index++;
			}
			
		}
		for (int i=0; i<this.lampdaRanges.size(); i++){
			int [][] count4Sel = tempHist[i];
			histT4SelAngleRange.add(count4Sel);
			
			// calculate average values
			averPhi[i] = averPhi[i]/foundNodes[i];
			averLampda[i] = averLampda[i]/foundNodes[i];
		}
		this.minMaxAverT4SelAngleRange.add(minLampda);
		this.minMaxAverT4SelAngleRange.add(minPhi);
		this.minMaxAverT4SelAngleRange.add(maxLampda);
		this.minMaxAverT4SelAngleRange.add(maxPhi);
		this.minMaxAverT4SelAngleRange.add(averLampda);
		this.minMaxAverT4SelAngleRange.add(averPhi);
		tMin =0;
//		switch (jSSTypeID){
//		case 0: // 'H'
////			col = (Color.magenta); break;
//		case 1: // 'S'
////			col = (Color.yellow); break;
//		case 2: // 'O'
////			col = (Color.cyan); break;
//		case 3: // 'A'
//		}
	}
	
	public void calcHistOfSphoxelMap(){
		ColorScaleView colView = this.contactView.getColorScaleView();
		HistogramView histView = this.contactView.getHistogramView();
		double [] scale;
		if (colView != null)
			scale = colView.getScaleValues();
		else 
			scale = histView.getScaleValues();
		this.histWholeAngleRange = new int[scale.length];
		this.hist4SelAngleRange = new Vector<int[]>();
		double phiRad, lampdaRad;
		double tMin, tMax, pMin, pMax;
		double ratio;
		int[][] counts = new int[this.lampdaRanges.size()][scale.length];
		for(int i=0;i<ratios.length;i++){
			for(int j=0;j<ratios[i].length;j++){	
				Iterator<Pair<Double>> itrP = this.lampdaRanges.iterator();
				Iterator<Pair<Double>> itrT = this.phiRanges.iterator();
				phiRad = (i*deltaRad);
				lampdaRad = (j*deltaRad);
				ratio = ratios[i][j];
				for(int k=0; k<scale.length-1; k++){
//					if (scale[k]==0)
//						k++;
					if (ratio>=scale[k] && ratio<scale[k+1]){
						histWholeAngleRange[k] += 1; // = histWholeAngleRange[k]+1;
						k = scale.length;
					}
				}
				int index = 0;
				while (itrP.hasNext()){
					Pair<Double> lampda = (Pair<Double>) itrP.next();
					Pair<Double> phi = (Pair<Double>) itrT.next();
					tMin = phi.getFirst();
					tMax = phi.getSecond();
					pMin = lampda.getFirst();
					pMax = lampda.getSecond();
					for(int k=0; k<scale.length-1; k++){
						if (phiRad>=tMin && phiRad<tMax && lampdaRad>=pMin && lampdaRad<pMax && ratio>=scale[k] && ratio<scale[k+1]){
							counts[index][k] += 1; //= counts[index][k]+1;
							k = scale.length;
						}
					}
					index++;
				}
			}
		}
		// copy values from int array to vector
		for (int i=0; i<this.lampdaRanges.size(); i++){
			int [] count4Sel = counts[i];
			hist4SelAngleRange.add(count4Sel);
		}
		
		ratio = 0;
	}
	
	/*------------------------ writing methods --------------------*/
	public void writeSphoxels(String filename){
		this.sphoxel.writeSphoxelOutput(filename);
	}
	public void writeTraces(String filename){
		this.nbhsTraces.writeNbhsTracesOutput(filename);
	}
	
	/*------------------------ drawing methods --------------------*/
	
	public GeneralPath rhombusShape(double xC, double yC, double width, double height){
		GeneralPath shape = new GeneralPath();
		shape.moveTo(xC-(width/2), yC);
		shape.lineTo(xC, yC-(height/2));
		shape.lineTo(xC+(width/2), yC);
		shape.lineTo(xC, yC+(height/2));
		shape.lineTo(xC-(width/2), yC);
		return shape;
	}
	
	public GeneralPath trapezium(Pair<Double> p1, Pair<Double> p2, Pair<Double> p3, Pair<Double> p4){
		GeneralPath shape = new GeneralPath();
		shape.moveTo(p1.getFirst(), p1.getSecond());
		shape.lineTo(p2.getFirst(), p2.getSecond());
		shape.lineTo(p3.getFirst(), p3.getSecond());
		shape.lineTo(p4.getFirst(), p4.getSecond());
		shape.lineTo(p1.getFirst(), p1.getSecond());
		return shape;
	}
	
	public GeneralPath trapezium(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
		GeneralPath shape = new GeneralPath();
		shape.moveTo(x1, y1);
		shape.lineTo(x2, y2);
		shape.lineTo(x3, y3);
		shape.lineTo(x4, y4);
		shape.lineTo(x1, y1);
		return shape;
	}
	
	public GeneralPath nGon(double[] x, double[] y){
		GeneralPath shape = new GeneralPath();
		if (x.length>0 && x.length==y.length){
			shape.moveTo(x[0], y[0]);
			for (int i=1; i<x.length; i++){
				shape.lineTo(x[i], y[i]);
			}
			shape.lineTo(x[0], y[0]);
		}
		return shape;
	}
	
	public GeneralPath pathLine(double[] x, double[] y){
		GeneralPath shape = new GeneralPath();
		if (x.length>0 && x.length==y.length){
			shape.moveTo(x[0], y[0]);
			for (int i=1; i<x.length; i++){
				shape.lineTo(x[i], y[i]);
			}
		}
		return shape;
	}
	
	public double getOffSet(double lampdaRad, double phiRad){
//		double lampdaRad = 0;
		double xPos=0;	
		if (this.mapProjType==kavrayskiyMapProj){
//			lampdaRad = translateXCoordRespective2Orig(lampdaRad);
			lampdaRad -= Math.PI;
			phiRad -= (Math.PI/2);
			
			xPos = lampda2Kavrayskiy(lampdaRad, phiRad);		
			xPos = ( (Math.PI+xPos) *this.voxelsize )+this.border;
		}
		else if (this.mapProjType==cylindricalMapProj){
			lampdaRad -= Math.PI;
			xPos = ( (Math.PI+lampdaRad) *this.voxelsize )+this.border;
		}
		return xPos;
	}
	
	private double distPointsOnSphere(double l1, double p1, double l2, double p2){
		double dist = 0;
//		double R = this.g2dSize.getHeight() / 2;
		
		double x1 = this.rSphere*Math.sin(p1)*Math.cos(l1);
		double y1 = this.rSphere*Math.sin(p1)*Math.sin(l1);
		double z1 = this.rSphere*Math.cos(p1);
		double x2 = this.rSphere*Math.sin(p2)*Math.cos(l2);
		double y2 = this.rSphere*Math.sin(p2)*Math.sin(l2);
		double z2 = this.rSphere*Math.cos(p2);
		
		dist = Math.sqrt(Math.pow(x1-x2,2) + Math.pow(y1-y2,2) + Math.pow(z1-z2,2));
		return dist;
	}
	
	private void drawSphoxelMap(Graphics2D g2d){
		double val;
//		double phiDeg, lampdaDeg; 
//		double phiRad, lampdaRad;
		Color col = Color.white; // = new Color(24,116,205,255);
		ColorScale scale = new ColorScale();
		//		double deltaDeg = 180/this.numSteps;
				
		// ----- color representation -----
		for(int i=0;i<ratios.length;i++){
			for(int j=0;j<ratios[i].length;j++){						
				val = ratios[i][j]; // add some scaling and shifting --> range 0:255				
				// ----- remove outliers
				double minRatio2Use = this.minRatio;
				double maxRatio2Use = this.maxRatio;
				if (this.removeOutliers){
					if (minRatio2Use<this.minAllowedRat)
						minRatio2Use = this.minAllowedRat;
					if (maxRatio2Use>this.maxAllowedRat)
						maxRatio2Use = this.maxAllowedRat;
					if (val<this.minAllowedRat)
						val = this.minAllowedRat;
					else if (val>this.maxAllowedRat)
						val = this.maxAllowedRat;
				}						
				// ----- compute alpha and set color	
				float alpha = (float) (Math.abs(val)/Math.abs(maxRatio2Use)); // if val>0
				if(val<0){
					alpha = (float) (Math.abs(val)/Math.abs(minRatio2Use));
					val = -val/minRatio2Use;
				}
				else
					val = val/maxRatio2Use;
				switch (this.chosenColourScale){
				case ContactStatusBar.BLUERED:
					col = scale.getColor4BlueRedScale(val,alpha); break;
				case ContactStatusBar.HOTCOLD:
					col = scale.getColor4HotColdScale(val,alpha); break;
				case ContactStatusBar.RGB:
					col = scale.getColor4RGBscalePolar((float)val, alpha, -1); break;
				}	
				
				g2d.setColor(col);	
			
				drawSphoxel(g2d, i, j);
			}
		}
	}
	
	private void drawSphoxel(Graphics2D g2d, int i, int j){
		Shape shape = null;
		double xPos1, yPos1, xPos2, yPos2, xPos3, yPos3, xPos4, yPos4;
		double phiRad = (i*deltaRad); // 0<phi<2*PI
		double lampdaRad = (j*deltaRad); // 0<lampda<2*PI
		// ---- create and draw shape (rectangle or trapezoid)		
		if (this.mapProjType==kavrayskiyMapProj){
			
			xPos1 = getScreenPosFromRad(lampdaRad, phiRad).getFirst();
			yPos1 = getScreenPosFromRad(lampdaRad, phiRad).getSecond();
			xPos2 = getScreenPosFromRad(lampdaRad+deltaRad, phiRad).getFirst();
			yPos2 = getScreenPosFromRad(lampdaRad+deltaRad, phiRad).getSecond();	
			phiRad += deltaRad;
			xPos3 = getScreenPosFromRad(lampdaRad+deltaRad, phiRad).getFirst();
			yPos3 = getScreenPosFromRad(lampdaRad+deltaRad, phiRad).getSecond();
			xPos4 = getScreenPosFromRad(lampdaRad, phiRad).getFirst();
			yPos4 = getScreenPosFromRad(lampdaRad, phiRad).getSecond();	
			phiRad -= deltaRad;
			
			if (xPos1>xPos2){
				double xPosUL=0, xPosLL=0, xPosUR=getSize().getWidth(), xPosLR=getSize().getWidth();
				xPosUL = getScreenPosFromRad(Math.PI-this.origCoordinates.getFirst(), phiRad).getFirst();
				xPosLL = getScreenPosFromRad(Math.PI-this.origCoordinates.getFirst(), phiRad+deltaRad).getFirst();
				xPosUR = getScreenPosFromRad(Math.PI-this.origCoordinates.getFirst()+2*Math.PI, phiRad).getFirst();
				xPosLR = getScreenPosFromRad(Math.PI-this.origCoordinates.getFirst()+2*Math.PI, phiRad+deltaRad).getFirst();
				shape = trapezium(xPos1, yPos1, xPosUR, yPos2, xPosLR, yPos3, xPos4, yPos4);	
				g2d.draw(shape);	
				g2d.fill(shape);
				shape = trapezium(xPosUL, yPos1, xPos2, yPos2, xPos3, yPos3, xPosLL, yPos4);	
				g2d.draw(shape);	
				g2d.fill(shape);
			}
			else 
			{
				shape = trapezium(xPos1, yPos1, xPos2, yPos2, xPos3, yPos3, xPos4, yPos4);	
				g2d.draw(shape);
				g2d.fill(shape);		
			}					
		}
		else if (this.mapProjType==cylindricalMapProj){				
			xPos1 = j*pixelWidth;
			yPos1 = i*pixelWidth;
			xPos1 = translateXPixelCoordRespective2Orig(xPos1) +this.border;
			yPos1 = translateYPixelCoordRespective2Orig(yPos1) +this.border;	
			shape = new Rectangle2D.Double(xPos1, yPos1, pixelWidth, pixelHeight);	
			g2d.draw(shape);
			g2d.fill(shape);
			
			if (xPos1+pixelWidth > this.xDim)
				xPos1 -= xDim;
			if (yPos1+pixelHeight > this.yDim)
				yPos1 -= yDim;
			shape = new Rectangle2D.Double(xPos1, yPos1, pixelWidth, pixelHeight);
			g2d.draw(shape);
			g2d.fill(shape);
		}
		else if (this.mapProjType==azimuthalMapProj){	
			xPos1 = getScreenPosFromRad(lampdaRad, phiRad).getFirst();
			yPos1 = getScreenPosFromRad(lampdaRad, phiRad).getSecond();	
			xPos2 = getScreenPosFromRad(lampdaRad+deltaRad, phiRad).getFirst();
			yPos2 = getScreenPosFromRad(lampdaRad+deltaRad, phiRad).getSecond();	
			xPos3 = getScreenPosFromRad(lampdaRad+deltaRad, phiRad+deltaRad).getFirst();
			yPos3 = getScreenPosFromRad(lampdaRad+deltaRad, phiRad+deltaRad).getSecond();	
			xPos4 = getScreenPosFromRad(lampdaRad, phiRad+deltaRad).getFirst();
			yPos4 = getScreenPosFromRad(lampdaRad, phiRad+deltaRad).getSecond();
			if (isOnFrontView(lampdaRad, phiRad))
				shape = trapezium(xPos1, yPos1, xPos2, yPos2, xPos3, yPos3, xPos4, yPos4);
			else
				shape = trapezium(xPos1+(2*this.rSphere), yPos1, xPos2+(2*this.rSphere), yPos2, xPos3+(2*this.rSphere), yPos3, xPos4+(2*this.rSphere), yPos4);		
			g2d.draw(shape);
			g2d.fill(shape);	
		}
	}
	
	@SuppressWarnings("unused")
	private void drawSphoxels(Graphics2D g2d){
		Shape shape = null;
		double xPos;
		double yPos;
		double val;
		Color col = Color.white; // = new Color(24,116,205,255);
		ColorScale scale = new ColorScale();
		
		g2d.setBackground(backgroundColor);
		shape = new Rectangle2D.Float(0, 0, this.pixelWidth*this.ratios[0].length, this.pixelHeight*this.ratios.length);
		g2d.setColor(backgroundColor);
		g2d.draw(shape);
		g2d.fill(shape);
				
		// ----- color representation -----
		for(int i=0;i<ratios.length;i++){
			for(int j=0;j<ratios[i].length;j++){
//				xPos = j*pixelWidth +this.border;
//				yPos = i*pixelHeight +this.border;
				xPos = translateXPixelCoordRespective2Orig(j*pixelWidth) +this.border;
				yPos = translateYPixelCoordRespective2Orig(i*pixelHeight) +this.border;				
//				shape = new Rectangle2D.Float(xPos+this.border, yPos+this.border+this.yBorderThres, pixelWidth, pixelHeight);
				shape = new Rectangle2D.Double(xPos, yPos, pixelWidth, pixelHeight);	
//				shape = trapezium(xPos, yPos, xPos+pixelWidth, yPos, xPos+pixelWidth, yPos+pixelHeight, xPos, yPos+pixelHeight);
				
				val = ratios[i][j]; // add some scaling and shifting --> range 0:255
				
				// ----- remove outliers
				double minRatio2Use = this.minRatio;
				double maxRatio2Use = this.maxRatio;
				if (this.removeOutliers){
//					if (this.minRatio<this.minAllowedRat)
//						this.minRatio = this.minAllowedRat;
//					if (this.maxRatio>this.maxAllowedRat)
//						this.maxRatio = this.maxAllowedRat;
					if (minRatio2Use<this.minAllowedRat)
						minRatio2Use = this.minAllowedRat;
					if (maxRatio2Use>this.maxAllowedRat)
						maxRatio2Use = this.maxAllowedRat;
					if (val<this.minAllowedRat)
						val = this.minAllowedRat;
					else if (val>this.maxAllowedRat)
						val = this.maxAllowedRat;
				}
								
				// ----- compute alpha and set color	
				float alpha = (float) (Math.abs(val)/Math.abs(maxRatio2Use)); // if val>0
				if(val<0){
					alpha = (float) (Math.abs(val)/Math.abs(minRatio2Use));
					val = -val/minRatio2Use;
				}
				else
					val = val/maxRatio2Use;
				if (alpha>1.0f) alpha=1.0f;
				if (alpha<0.0f) alpha=0.0f;
//				val = -val;
				switch (this.chosenColourScale){
				case ContactStatusBar.BLUERED:
					col = scale.getColor4BlueRedScale(val,alpha); break;
				case ContactStatusBar.HOTCOLD:
					col = scale.getColor4HotColdScale(val, alpha);
				case ContactStatusBar.RGB:
					col = scale.getColor4RGBscalePolar((float)val, alpha, -1);
				}
				
				g2d.setColor(col);
//				g2d.setColor(Color.black); // Test purpose
				g2d.draw(shape);
//				g2d.setColor(col);
				g2d.fill(shape);
				if (xPos+pixelWidth > this.xDim)
					xPos -= xDim;
				if (yPos+pixelHeight > this.yDim)
					yPos -= yDim;
				shape = new Rectangle2D.Double(xPos, yPos, pixelWidth, pixelHeight);
				g2d.draw(shape);
//				g2d.setColor(col);
				g2d.fill(shape);
			}
		}
	}
	
	public Color getNodeColor4SSType(int jSSTypeID){
		Color col = Color.black;
		switch (jSSTypeID){
		case 0: // 'H'
			col = (helixSSTColor); break;
		case 1: // 'S'
			col = (sheetSSTColor); break;
		case 2: // 'O'
			col = (otherSSTColor); break;
		}
		return col;
	}
	
	private void drawNBHSNode(Graphics2D g2d, boolean specialRes, float[] node){
		GeneralPath rhombus = null;
		Shape circle = null;
		int iNum, jNum, jResID, jSSType;
		double xPos, yPos;
		double phiRad, lampdaRad;
		Font f = new Font("Dialog", Font.PLAIN, 12);
		float radius = 3.f;
		String nodeName;
		
		iNum = (int) node[1];
		jNum = (int) node[2];
		phiRad = node[3];
		lampdaRad = node[4];
		jResID = (int) node[5];
		jSSType = (int) node[6];
		
		lampdaRad += Math.PI;
		
		xPos = getScreenPosFromRad(lampdaRad, phiRad).getFirst();
		yPos = getScreenPosFromRad(lampdaRad, phiRad).getSecond();
		if (this.mapProjType==azimuthalMapProj && !isOnFrontView(lampdaRad, phiRad))
			xPos += 2*this.rSphere;		
		
		// ---- draw geometric object for each residue
		Color col = getNodeColor4SSType(jSSType);
		g2d.setColor(col);
		if (specialRes){
			radius = 6.f;
			// create rhombus for residue types contained in nbhstring
			rhombus = rhombusShape(xPos, yPos+this.yBorderThres, 2*radius, 2*radius);
			g2d.draw(rhombus);
			g2d.fill(rhombus);
			f = new Font("Dialog", Font.PLAIN, 14);
		}
		else {
			radius = 3.f;
			// create ellipse for residue
			circle = new Ellipse2D.Double( xPos-radius, yPos-radius+this.yBorderThres,2*radius, 2*radius);
			g2d.draw(circle);
			g2d.fill(circle);
			f = new Font("Dialog", Font.PLAIN, 12);
		}
		g2d.setFont(f);
		nodeName = Character.toString(aas[jResID]) +" "+ String.valueOf(jNum-iNum);
		if (specialRes)
			g2d.setColor(Color.black);
		else
			g2d.setColor(new Color(70,70,70));
		g2d.drawString(nodeName, (float)(xPos+radius), (float)(yPos+radius+this.yBorderThres));
	}
	
	private void drawNBHSEdges(Graphics2D g2d, float[] node, int nodeID, int lineID, int j){
//		Shape line = null;
		int gID, iNum, jNum, gIdNB, iNumNB, jNumNB; 
		float[] nbNode;
		ColorScale scale = new ColorScale();
		Color col = null;		
		double phiRad, lampdaRad, phiRadNB, lampdaRadNB;
		
		gID = (int) node[0];
		iNum = (int) node[1];
		jNum = (int) node[2];
		phiRad = node[3];
		lampdaRad = node[4];				
		lampdaRad += Math.PI; // cause lampda is expected to have a value of [0:pi]
		
		// --- gradient color edges between connected nodes
		nbNode = (float[]) this.nbhsNodes.get(j);
		phiRadNB = nbNode[3];
		lampdaRadNB = nbNode[4];
		lampdaRadNB += Math.PI; // cause lampda is expected to have a value of [0:pi]
		
		if (node[0]==nbNode[0] && node[1]==nbNode[1]){
//			System.out.print(nodeID +"\t" + lineID +"\t" + this.numNodesPerLine[lineID] +"\t");
			float ratio = (float)(nodeID+1)/(float)this.numNodesPerLine[lineID];
//			ratio = ((this.maxDistsJRes[lineID]-this.minDistsJRes[lineID])*node[2] + this.minDistsJRes[lineID])/this.maxDistsJRes[lineID];
//			System.out.print(ratio +"\t");
			
			// --- use ratio for color scale
//			g2d.setColor(scale.getColor4YellowRedScale(ratio, 1.0f));
			g2d.setColor(scale.getColor4RGBscale(ratio, 1.0f, 1));
			
			gIdNB = (int) nbNode[0];
			iNumNB = (int) nbNode[1];
			jNumNB = (int) nbNode[2];
			if (gID==gIdNB && iNum==iNumNB){
//				// -- equal scaling
//				if (jNum-iNum < 0)
//					ratio = -1 + (float)(jNum-this.minDistsJRes[lineID])/(float)(iNum-this.minDistsJRes[lineID]);
//				else 
//					ratio = (float)(jNum-iNum)/(float)(this.maxDistsJRes[lineID]-iNum);
				
//				// -- logarithmic scaling
//				if (jNum-iNum < 0)
//					ratio = -1 + (float)Math.log(jNum-this.minDistsJRes[lineID]+1)/(float)Math.log(iNum-this.minDistsJRes[lineID]+1);
////					ratio = -1 + (float)Math.log10(jNum-this.minDistsJRes[lineID]+1)/(float)Math.log10(iNum-this.minDistsJRes[lineID]+1);
//				else 
//					ratio = (float)Math.log(jNum-iNum+1)/(float)Math.log(this.maxDistsJRes[lineID]-iNum+1);
////					ratio = (float)Math.log10(jNum-iNum+1)/(float)Math.log10(this.maxDistsJRes[lineID]-iNum+1);
//				System.out.print(iNum+"_"+jNum+"_"+this.minDistsJRes[lineID]+"_"+this.maxDistsJRes[lineID]+":"+ratio +"\t");
//				g2d.setColor(scale.getColor4RGBscalePolar(ratio, 1.0f, 1));
				
				// -- shortrange scaling: |jNum-iNum|>ShortRangeThreshold --> blue
				int thres1 = 9; // 1-9:short range  9-25:middle range  25-n/9-n:long range
				int thres2 = 25;
				col = Color.black;
				if (Math.abs(jNum-iNum)<=thres1){
					ratio = +1 * (float)Math.abs(jNum-iNum)/(float)(thres1);
					// scale on range 0.2:0.8
					ratio = 0.2f + (ratio*(0.8f-0.2f));
					col = scale.getColor4GreyValueRange(ratio, 1);
				}
				else if (Math.abs(jNum-iNum)<=thres2){
					if (jNum-iNum < 0)
						ratio = -1 * (float)(Math.abs(jNum-iNum)-thres1)/(float)(thres2-thres1);
					else 
						ratio = +1 * (float)(Math.abs(jNum-iNum)-thres1)/(float)(thres2-thres1);
					col = scale.getColor4HotColdScale(ratio, 1.0f);
				}
				else {
					if (jNum-iNum < 0)
						ratio = -1.0f;
					else 
						ratio = +1.0f;
					col = scale.getColor4HotColdScale(ratio, 1.0f);
				}
//				System.out.print(iNum+"_"+jNum+"_"+Math.abs(jNum-iNum)+":"+ratio +"\t");
				g2d.setColor(col);					
				
				
				if (iNum>jNum && iNum<jNumNB && this.paintCentralResidue){
					// paint central residue
//					circle = new Ellipse2D.Float( xPosINum-radius, yPosINum+this.yBorderThres-radius,2*radius, 2*radius);
					
					drawNBHSEdge(g2d, lampdaRad, phiRad, Math.PI, Math.PI/2);
					drawNBHSEdge(g2d, Math.PI, Math.PI/2, lampdaRadNB, phiRadNB);
				}
				else {
					drawNBHSEdge(g2d, lampdaRad, phiRad, lampdaRadNB, phiRadNB);
				}	
			}									
			nodeID++;
		}		
	}
	
	private void drawNBHSEdge(Graphics2D g2d, double lampdaRad, double phiRad, double lampdaRadNB, double phiRadNB){
		Shape line = null;
		double xPos, yPos, xPosNB=0, yPosNB=0;
		xPos = getScreenPosFromRad(lampdaRad, phiRad).getFirst();
		yPos = getScreenPosFromRad(lampdaRad, phiRad).getSecond();
		if (this.mapProjType==azimuthalMapProj && !isOnFrontView(lampdaRad, phiRad))
			xPos += (2*this.rSphere);
		xPosNB = getScreenPosFromRad(lampdaRadNB, phiRadNB).getFirst();
		yPosNB = getScreenPosFromRad(lampdaRadNB, phiRadNB).getSecond();
		if (this.mapProjType==azimuthalMapProj && !isOnFrontView(lampdaRadNB, phiRadNB))
			xPosNB += (2*this.rSphere);
		
		// -- test if edge should be wrapped
		if (this.mapProjType==cylindricalMapProj && Math.abs(xPosNB-xPos)>this.xDim/2){
			if (xPos<xPosNB){
				line = new Line2D.Double(xPos, yPos+this.yBorderThres, xPosNB-this.xDim, yPosNB+this.yBorderThres);		
				g2d.draw(line);
				line = new Line2D.Double(xPos+this.xDim, yPos+this.yBorderThres, xPosNB, yPosNB+this.yBorderThres);		
				g2d.draw(line);
			}
			else{
				line = new Line2D.Double(xPos-this.xDim, yPos+this.yBorderThres, xPosNB, yPosNB+this.yBorderThres);		
				g2d.draw(line);
				line = new Line2D.Double(xPos, yPos+this.yBorderThres, xPosNB+this.xDim, yPosNB+this.yBorderThres);		
				g2d.draw(line);
			}
		}
		else if (this.mapProjType==kavrayskiyMapProj && Math.abs(xPosNB-xPos)>this.xDim/4){
			double phiRadM = phiRad + (phiRadNB-phiRad)/2;
//			double test = this.origCoordinates.getFirst()-Math.PI;
			double borderXPos = getScreenPosFromRad(Math.PI-this.origCoordinates.getFirst(), phiRadM).getFirst();
			double borderYPos = getScreenPosFromRad(lampdaRadNB, phiRadM).getSecond();
			if (xPos<xPosNB){
				line = new Line2D.Double(xPos, yPos+this.yBorderThres, borderXPos, borderYPos+this.yBorderThres);		
				g2d.draw(line);
				borderXPos = getScreenPosFromRad((3*Math.PI)-this.origCoordinates.getFirst(), phiRadM).getFirst();
				line = new Line2D.Double(xPosNB, yPosNB+this.yBorderThres, borderXPos, borderYPos+this.yBorderThres);		
				g2d.draw(line);							
			}
			else {
				line = new Line2D.Double(xPosNB, yPosNB+this.yBorderThres, borderXPos, borderYPos+this.yBorderThres);		
				g2d.draw(line);	
				borderXPos = getScreenPosFromRad((3*Math.PI)-this.origCoordinates.getFirst(), phiRadM).getFirst();	
				line = new Line2D.Double(xPos, yPos+this.yBorderThres, borderXPos, borderYPos+this.yBorderThres);		
				g2d.draw(line);					
			}
		}
		else if (this.mapProjType==azimuthalMapProj){
			if (isOnFrontView(lampdaRad, phiRad) != isOnFrontView(lampdaRadNB, phiRadNB)){
				double borderXFront, borderYFront, borderXBack, borderYBack;
				// compute border positions of front and back view
				Pair<Double> intermP = intermediatePointAzim(lampdaRad, phiRad, lampdaRadNB, phiRadNB);
				borderXFront = getScreenPosFromRad(intermP.getFirst(), intermP.getSecond()).getFirst(); //intermP.getFirst();
				borderYFront = getScreenPosFromRad(intermP.getFirst(), intermP.getSecond()).getSecond(); //intermP.getSecond();
				borderXBack = borderXFront + (2*this.rSphere);
				borderYBack = borderYFront;
				if (isOnFrontView(lampdaRad, phiRad)){
					line = new Line2D.Double(xPos, yPos+this.yBorderThres, borderXFront, borderYFront+this.yBorderThres);		
					g2d.draw(line);
					line = new Line2D.Double(xPosNB, yPosNB+this.yBorderThres, borderXBack, borderYBack+this.yBorderThres);		
					g2d.draw(line);
				}
				else {
					line = new Line2D.Double(xPos, yPos+this.yBorderThres, borderXBack, borderYBack+this.yBorderThres);		
					g2d.draw(line);
					line = new Line2D.Double(xPosNB, yPosNB+this.yBorderThres, borderXFront, borderYFront+this.yBorderThres);		
					g2d.draw(line);
				}				
			}
			else 
			{
				// standard
				line = new Line2D.Double(xPos, yPos+this.yBorderThres, xPosNB, yPosNB+this.yBorderThres);		
				g2d.draw(line);
			}
		}
		else 
		{
			line = new Line2D.Double(xPos, yPos+this.yBorderThres, xPosNB, yPosNB+this.yBorderThres);		
			g2d.draw(line);
		}
	}
	
	/* Computes the intermediate point top wrap line between two points, if points are situated on different views.
	 * The values of for lampda and phi are expected to be within ranges [0:2PI] and [0:PI]
	 * */
	private Pair<Double> intermediatePointAzim(double l1, double p1, double l2, double p2){
		Pair<Double> lp;		
		double l, p, lS=0, lE=0, dL = 0.1;
		if (isOnFrontView(l1, p1)){
			lS = l1; p = p1;
			lE = l2;
		}
		else{
			lS = l2; p = p2;
			lE  = l1;
		}
		if (lS>lE)
			dL = -1*dL;
		l=lS;
		lp = new Pair<Double>(l,p);
		while ( ((l>=lS && l<=lE) || (l>=lE && l<=lS)) && 
				distPointsOnSphere(this.centerOfProjection.getFirst(), this.centerOfProjection.getSecond(), l, p)<this.maxDistPoints ){
			lp = new Pair<Double>(l,p);
			l += dL;
			p = ((p2-p1)/(l2-l1))*(l-l1) + p1;
		}
		return lp;
	}
	
	private void drawNBHSTraces(Graphics2D g2d){
		int gID, iNum, jResID;
		float[] node, nbNode;
		int lineID = 0;
		int nodeID = 0;

		int nbhsIndexC = 0;
		boolean specialRes = false;
				
		for(int i=0; i<this.nbhsNodes.size(); i++){
			
			node = (float[]) this.nbhsNodes.get(i);
//			System.out.println("0:graphi_id + '\t' + 1:i_num + '\t' + 2:j_num + '\t' + 3:phi + '\t' + 4:lampda");
			gID = (int) node[0];
			iNum = (int) node[1];
			jResID = (int) node[5];
			
			// ---- check if jResType is element of nbhstring
			specialRes = false;
			if (i>0){
				nbNode = this.nbhsNodes.get(i-1);
				if (nbNode[0]!=gID || nbNode[1]!=iNum)
					nbhsIndexC = 0;				
			}
			// simple differentiation
			if (this.nbhString.contains(String.valueOf(aas[jResID]))){
				specialRes = true;
			}
			// ordering is of importance
			if (nbhsIndexC > 0 && aas[jResID]==this.nbhsRes[nbhsIndexC-1])
				specialRes = true;
			if (aas[jResID] == this.nbhsRes[nbhsIndexC]){
				specialRes = true;
				nbhsIndexC++;
			}				
			
			// ---- draw geometric object for each residue
			drawNBHSNode(g2d, specialRes, node);				
			
			// jump over empty traces
			while (this.numNodesPerLine[lineID]==0){
				lineID++;
				nodeID = 0;
			}
			int j=i+1;
			
			if (j<this.nbhsNodes.size()){				
				drawNBHSEdges(g2d, node, nodeID, lineID, j);
			}	
			else {
				lineID++;
//			    System.out.println();
				nodeID = 0;
			}		
		}	
		
	}
	
	private void drawLongitude(Graphics2D g2d, double lampdaRad, double phiRad){
		double xS, xE, yS, yE;
		xS = getScreenPosFromRad(lampdaRad, phiRad).getFirst();
		xE = getScreenPosFromRad(lampdaRad, phiRad+deltaRad).getFirst();
		yS = getScreenPosFromRad(lampdaRad, phiRad).getSecond();
		yE = getScreenPosFromRad(lampdaRad, phiRad+deltaRad).getSecond();
		if (this.mapProjType==azimuthalMapProj){
			if (isOnFrontView(lampdaRad, phiRad)==isOnFrontView(lampdaRad+deltaRad, phiRad)){
				if (!isOnFrontView(lampdaRad, phiRad)){
					xS += (2*this.rSphere);
					xE += (2*this.rSphere);
				}
				drawThickVerticalLine(g2d, xS, yS, xE, yE);
			}
		}
		else 
			drawThickVerticalLine(g2d, xS, yS, xE, yE);
	}
	
	private void drawLatitude(Graphics2D g2d, double lampdaRad, double phiRad){
		double xS, xE, yS, yE;
		xS = getScreenPosFromRad(lampdaRad, phiRad).getFirst();
		xE = getScreenPosFromRad(lampdaRad+deltaRad, phiRad).getFirst();
		yS = getScreenPosFromRad(lampdaRad, phiRad).getSecond();
		yE = getScreenPosFromRad(lampdaRad+deltaRad, phiRad).getSecond();
		if (this.mapProjType==azimuthalMapProj && (isOnFrontView(lampdaRad, phiRad)==isOnFrontView(lampdaRad+deltaRad, phiRad))){
			if (!isOnFrontView(lampdaRad, phiRad)){
				xS += (2*this.rSphere);
				xE += (2*this.rSphere);
			}
			drawThickHorizontalLine(g2d, xS, yS, xE, yE);
		}
//		else
//			drawThickHorizontalLine(g2d, xS, yS, xE, yE);
	}
	
	private void drawThickVerticalLine(Graphics2D g2d, double xS, double yS, double xE, double yE){
		Shape line;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
		xS=xS-dL; xE=xE-dL;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);	
	}
	
	private void drawThickHorizontalLine(Graphics2D g2d, double xS, double yS, double xE, double yE){
		Shape line;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
		yS = yS-dL; yE = yE-dL;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
	}
	
//	private void drawThickLine(Graphics2D g2d, double xS, double yS, double xE, double yE){
//		Shape line;
//		line = new Line2D.Double(xS,yS,xE,yE);
//		g2d.draw(line);
//		xS=xS-dL; xE=xE-dL;
//		yS = yS-dL; yE = yE-dL;
//		line = new Line2D.Double(xS,yS,xE,yE);
//		g2d.draw(line);
//	}
	
	private void drawLongitudes(Graphics2D g2d){
		// Laengengrad
		g2d.setColor(this.longitudeColor);
		double xS, xE, yS, yE;
		double phiRad, lampdaRad;
		
		
		if (this.mapProjType==kavrayskiyMapProj || this.mapProjType==azimuthalMapProj){			
			for(int i=0;i<ratios.length;i++){
				phiRad = (i*deltaRad);
				lampdaRad = 0.0;
				drawLongitude(g2d, lampdaRad, phiRad);					
				lampdaRad = Math.PI/2;
				drawLongitude(g2d, lampdaRad, phiRad);					
				lampdaRad = Math.PI;
				drawLongitude(g2d, lampdaRad, phiRad);
				lampdaRad = 3*Math.PI/2;
				drawLongitude(g2d, lampdaRad, phiRad);				
				lampdaRad = 2*Math.PI;
				drawLongitude(g2d, lampdaRad, phiRad);
			}
		}
		else if (this.mapProjType==cylindricalMapProj){
			yS = ((0.0) *this.voxelsize) +this.border;
			yE = (Math.PI *this.voxelsize) +this.border;
			
			xS = (translateXCoordRespective2Orig(0.0) *this.voxelsize) +this.border;  xE=xS;
			drawThickVerticalLine(g2d, xS, yS, xE, yE);
			xS = (translateXCoordRespective2Orig(Math.PI/2) *this.voxelsize) +this.border;  xE=xS;
			drawThickVerticalLine(g2d, xS, yS, xE, yE);
			xS = (translateXCoordRespective2Orig(Math.PI) *this.voxelsize) +this.border;  xE=xS;
			drawThickVerticalLine(g2d, xS, yS, xE, yE);
			xS = (translateXCoordRespective2Orig(3*Math.PI/2) *this.voxelsize) +this.border;  xE=xS;
			drawThickVerticalLine(g2d, xS, yS, xE, yE);
			xS = (translateXCoordRespective2Orig(2*Math.PI) *this.voxelsize) +this.border;  xE=xS;
			drawThickVerticalLine(g2d, xS, yS, xE, yE);
		}
//		else if (this.mapProjType==azimuthalMapProj){
//			
//		}
					
	}
	
	private void drawLatitudes(Graphics2D g2d){
		// Breitengrad
		g2d.setColor(this.latitudeColor);
		double xS, xE, yS, yE; //yE = yS;
		
		if (this.mapProjType==azimuthalMapProj){
			double phiRad, lampdaRad;			
			for(int j=0;j<ratios[0].length;j++){
				lampdaRad = (j*deltaRad);
				
				phiRad = Math.PI/4;
				drawLatitude(g2d, lampdaRad, phiRad);
				phiRad = Math.PI/2;
				drawLatitude(g2d, lampdaRad, phiRad);
				phiRad = 3*Math.PI/4;
				drawLatitude(g2d, lampdaRad, phiRad);
			}
		}
		else if (this.mapProjType==kavrayskiyMapProj || this.mapProjType==cylindricalMapProj){		
			xS = ((0.0f) *this.voxelsize) +this.border;
			xE = (((2*Math.PI)) *this.voxelsize) +this.border;
			yS = (translateYCoordRespective2Orig(Math.PI/2) *this.voxelsize) +this.border; yE = yS;
			drawThickHorizontalLine(g2d, xS, yS, xE, yE);
			yS = (translateYCoordRespective2Orig(Math.PI/4) *this.voxelsize) +this.border; yE = yS;
			drawThickHorizontalLine(g2d, xS, yS, xE, yE);
			yS = (translateYCoordRespective2Orig(3*Math.PI/4) *this.voxelsize) +this.border; yE = yS;
			drawThickHorizontalLine(g2d, xS, yS, xE, yE);		
		}
	}
	
	private void drawLongLatCentre(Graphics2D g2d){
		double bcenterx=0, bcentery=0;
		if (this.mapProjType==kavrayskiyMapProj){
			bcenterx = getScreenPosFromRad(Math.PI, Math.PI/2).getFirst();
			bcentery = getScreenPosFromRad(Math.PI, Math.PI/2).getSecond();
		}
		else if (this.mapProjType==cylindricalMapProj){
			bcenterx = (this.origCoordinates.getFirst() * this.voxelsize) + this.border;
			bcentery = (this.origCoordinates.getSecond() * this.voxelsize) + this.border;
		}
		else if (this.mapProjType==azimuthalMapProj){
			bcenterx = getScreenPosFromRad(Math.PI, Math.PI/2).getFirst();
			bcentery = getScreenPosFromRad(Math.PI, Math.PI/2).getSecond();
//			if (distPointsOnSphere(Math.PI, Math.PI/2, this.centerOfProjection.getFirst(), this.centerOfProjection.getSecond()) > this.maxDistPoints)
			if (!isOnFrontView(Math.PI, Math.PI/2))
				bcenterx += (2*this.rSphere);		
		}
		
		Shape circle = new Ellipse2D.Double(bcenterx-30, bcentery-30, 60, 60);
		g2d.draw(circle);
		circle = new Ellipse2D.Double(bcenterx-29, bcentery-29, 58, 58);
		g2d.draw(circle);		

//		bcenterx = getScreenPosFromRad(this.centerOfProjection.getFirst(), this.centerOfProjection.getSecond()).getFirst();
//		bcentery = getScreenPosFromRad(this.centerOfProjection.getFirst(), this.centerOfProjection.getSecond()).getSecond();
//		g2d.setColor(Color.blue);
//		circle = new Ellipse2D.Double(bcenterx-30, bcentery-30, 60, 60);
//		g2d.draw(circle);
	}
	
	private void drawCrosshair(Graphics2D g2d){
		// only in case of range selection we draw a diagonal cursor
		// drawing the cross-hair
		g2d.setColor(crosshairColor);
		g2d.drawLine(mousePos.x, 0, mousePos.x, g2dSize.height);
		g2d.drawLine(0, mousePos.y, g2dSize.width, mousePos.y);	
		g2d.drawLine(mousePos.x-1, 0, mousePos.x-1, g2dSize.height);
		g2d.drawLine(0, mousePos.y-1, g2dSize.width, mousePos.y-1);	
	}

//	private void drawRulerCrosshair(Graphics2D g2d) {
//		int x1,x2,y1,y2;
//		g2d.setColor(crosshairColor);
//		if(currentRulerMouseLocation == ResidueRuler.TOP || currentRulerMouseLocation == ResidueRuler.BOTTOM) {
//			x1 = currentRulerMousePos;
//			x2 = currentRulerMousePos;
//			y1 = 0;
//			y2 = outputSize;
//		} else {
//			x1 = 0;
//			x2 = outputSize;
//			y1 = currentRulerMousePos;
//			y2 = currentRulerMousePos;			
//		}
//		g2d.drawLine(x1, y1, x2, y2);
//	}
	
	private void drawOccupiedAngleRanges(Graphics2D g2d){
		Shape shape;
		double xPos1, yPos1, xPos2, yPos2, xPos3, yPos3, xPos4, yPos4;
		double xPosUL=0, xPosLL=0, xPosUR=getSize().getWidth(), xPosLR=getSize().getWidth();
		double yPosUL=0, yPosLL=0, yPosUR=getSize().getWidth(), yPosLR=getSize().getWidth();
		Iterator<Pair<Double>> itrP = this.lampdaRanges.iterator();
		Iterator<Pair<Double>> itrT = this.phiRanges.iterator();
		while (itrP.hasNext()){
			g2d.setColor(selAngleRangeColor);
			Pair<Double> lampda = (Pair<Double>) itrP.next();
			Pair<Double> phi = (Pair<Double>) itrT.next();
			
			double xS = lampda.getFirst();
			double xE = lampda.getSecond();
			double yS = phi.getFirst();
			double yE = phi.getSecond();
			xPos1 = getScreenPosFromRad(xS, yS).getFirst();
			yPos1 = getScreenPosFromRad(xS, yS).getSecond();
			xPos2 = getScreenPosFromRad(xE, yS).getFirst();
			yPos2 = getScreenPosFromRad(xE, yS).getSecond();
			xPos3 = getScreenPosFromRad(xE, yE).getFirst();
			yPos3 = getScreenPosFromRad(xE, yE).getSecond();
			xPos4 = getScreenPosFromRad(xS, yE).getFirst();
			yPos4 = getScreenPosFromRad(xS, yE).getSecond();
			if (this.mapProjType==azimuthalMapProj){
				if (!isOnFrontView(xS, yS))
					xPos1 += (2*this.rSphere);
				if (!isOnFrontView(xE, yS))
					xPos2 += (2*this.rSphere);
				if (!isOnFrontView(xE, yE))
					xPos3 += (2*this.rSphere);
				if (!isOnFrontView(xS, yE))
					xPos4 += (2*this.rSphere);
				
				if (isOnFrontView(xS, yS) != isOnFrontView(xE, yE)){
					// compute border positions of front and back view
					Pair<Double> intermP = intermediatePointAzim(xS, yS, xE, yS);
					xPosUL = getScreenPosFromRad(intermP.getFirst(), intermP.getSecond()).getFirst(); //intermP.getFirst();
					yPosUL = getScreenPosFromRad(intermP.getFirst(), intermP.getSecond()).getSecond(); //intermP.getSecond();
					xPosUR = xPosUL + (2*this.rSphere);
					yPosUR = yPosUL;
					intermP = intermediatePointAzim(xS, yE, xE, yE);
					xPosLL = getScreenPosFromRad(intermP.getFirst(), intermP.getSecond()).getFirst(); //intermP.getFirst();
					yPosLL = getScreenPosFromRad(intermP.getFirst(), intermP.getSecond()).getSecond(); //intermP.getSecond();
					xPosLR = xPosLL + (2*this.rSphere);
					yPosLR = yPosLL;
					double[] x;
					double[] y; 
					if (isOnFrontView(xS, yS)){
						x = new double[]{xPosUL, xPos1+1, xPos4, xPosLL};
						y = new double[]{yPosUL, yPos1, yPos4, yPosLL}; 
						shape = pathLine(x, y);
						g2d.draw(shape);
						x = new double[]{xPosUR, xPos2, xPos3, xPosLR};
						y = new double[]{yPosUR, yPos2, yPos3, yPosLR}; 
						shape = pathLine(x, y);
						g2d.draw(shape);
					}
					else {
						x = new double[]{xPosUL, xPos2, xPos3, xPosLL};
						y = new double[]{yPosUL, yPos2, yPos3, yPosLL}; 
						shape = pathLine(x, y);
						g2d.draw(shape);
						x = new double[]{xPosUR, xPos1, xPos4, xPosLR};
						y = new double[]{yPosUR, yPos1, yPos4, yPosLR}; 
						shape = pathLine(x, y);
						g2d.draw(shape);
					}				
				}
				
			}	
			if (this.mapProjType!=azimuthalMapProj && xPos1>xPos2){
//				System.out.println("Umklappen!!!!");
				xPosUL = getScreenPosFromRad(Math.PI-this.origCoordinates.getFirst(), yS).getFirst();
				xPosLL = getScreenPosFromRad(Math.PI-this.origCoordinates.getFirst(), yE).getFirst();
				xPosUR = getScreenPosFromRad(Math.PI-this.origCoordinates.getFirst()+2*Math.PI, yS).getFirst();
				xPosLR = getScreenPosFromRad(Math.PI-this.origCoordinates.getFirst()+2*Math.PI, yE).getFirst();
				double[] x = new double[]{xPosUR+1, xPos1, xPos4, xPosLR+1};
				double[] y = new double[]{yPos2, yPos1, yPos4, yPos3}; 
				shape = pathLine(x, y);
				g2d.draw(shape);
				x = new double[]{xPosUR+1, xPos1+1, xPos4+1, xPosLR+1};
				y = new double[]{yPos2+1, yPos1+1, yPos4-1, yPos3-1}; 
				shape = pathLine(x, y);
				g2d.draw(shape);
				x = new double[]{xPosUL-1, xPos2, xPos3, xPosLL-1};
				y = new double[]{yPos1, yPos2, yPos3, yPos4}; 
				shape = pathLine(x, y);
				g2d.draw(shape);
				x = new double[]{xPosUL-1, xPos2+1, xPos3+1, xPosLL-1};
				y = new double[]{yPos1+1, yPos2+1, yPos3-1, yPos4-1}; 
				shape = pathLine(x, y);
				g2d.draw(shape);				
			}
			else{
				shape = trapezium(xPos1, yPos1, xPos2, yPos2, xPos3, yPos3, xPos4, yPos4);			
				g2d.draw(shape);
				shape = trapezium(xPos1+1, yPos1+1, xPos2-1, yPos2+1, xPos3-1, yPos3-1, xPos4+1, yPos4-1);		
				g2d.draw(shape);				
			}
		}
		
	}
	
	private void drawSelectedAngleRange(Graphics2D g2d){
		if (dragging && contactView.getGUIState().getSelectionMode()==ContactGUIState.SelMode.RECT) {
			double xS, yS, xE, yE;
			double xPos1, yPos1, xPos2, yPos2, xPos3, yPos3, xPos4, yPos4;
			g2d.setColor(squareSelColor);
			
			Pair<Double> posP = screen2spherCoord(mousePressedPos);
			Pair<Double> posD = screen2spherCoord(mouseDraggingPos);
			xS = posP.getFirst();
			yS = posP.getSecond();
			xE = posD.getFirst();
			yE = posD.getSecond();
			if (this.mapProjType==kavrayskiyMapProj){
				xS = lampdaFromKavrayskiy(xS-Math.PI, yS-(Math.PI/2)) + Math.PI;
				xE = lampdaFromKavrayskiy(xE-Math.PI, yE-(Math.PI/2)) + Math.PI;
			}
			xS -= (this.origCoordinates.getFirst()-Math.PI);
			xE -= (this.origCoordinates.getFirst()-Math.PI);
			
			xPos1 = getScreenPosFromRad(xS, yS).getFirst();
			yPos1 = getScreenPosFromRad(xS, yS).getSecond();
			xPos2 = getScreenPosFromRad(xE, yS).getFirst();
			yPos2 = getScreenPosFromRad(xE, yS).getSecond();
			xPos3 = getScreenPosFromRad(xE, yE).getFirst();
			yPos3 = getScreenPosFromRad(xE, yE).getSecond();
			xPos4 = getScreenPosFromRad(xS, yE).getFirst();
			yPos4 = getScreenPosFromRad(xS, yE).getSecond();	
			if (this.mapProjType==azimuthalMapProj){
				if (!isOnFrontView(xS, yS))
					xPos1 += (2*this.rSphere);
				if (!isOnFrontView(xE, yS))
					xPos2 += (2*this.rSphere);
				if (!isOnFrontView(xE, yE))
					xPos3 += (2*this.rSphere);
				if (!isOnFrontView(xS, yE))
					xPos4 += (2*this.rSphere);
			}	
			Shape shape = trapezium(xPos1, yPos1, xPos2, yPos2, xPos3, yPos3, xPos4, yPos4);			
			g2d.draw(shape);
		} 
	}
	
	private void updateAnglePanel(){
		this.contStatBar.getAnglePanel().setIRes(String.valueOf(this.iRes));
		this.contStatBar.getAnglePanel().setJRes(String.valueOf(this.jRes));
		this.contStatBar.getAnglePanel().setISSType(String.valueOf(this.iSSType));
		this.contStatBar.getAnglePanel().setJSSType(String.valueOf(this.jSSType));
		this.contStatBar.getAnglePanel().setNBHS(this.nbhString);
		
		float val = (float)(Math.round(this.tmpLampdaRange.getFirst()*100))/100;
		this.contStatBar.getAnglePanel().setLampdaMin(String.valueOf(val));	
		val = (float)(Math.round(this.tmpAngleComb.getFirst()*100))/100;
		this.contStatBar.getAnglePanel().setLampdaMax(String.valueOf(val));	
		val = (float)(Math.round(this.tmpPhiRange.getFirst()*100))/100;
		this.contStatBar.getAnglePanel().setPhiMin(String.valueOf(val));	
		val = (float)(Math.round(this.tmpAngleComb.getSecond()*100))/100;
		this.contStatBar.getAnglePanel().setPhiMax(String.valueOf(val));
		
		this.contStatBar.repaint();		
	}
	
	// end drawing methods
	
	// --- Methods for coordinate handling (spherical and cartesian coordinates of projection) -----
	
	public boolean angleWithinValidRange(double xPos, double yPos){
		boolean valid = true;
		if (this.mapProjType==kavrayskiyMapProj){
			if (xPos<0 || xPos>2*Math.PI)
				valid = false;
		}
		else if (this.mapProjType==azimuthalMapProj){
			
		}
		return valid;
	}
	
	/*
	 * Pseudocylindrical projection of angle values
	 * Expects values for phi[-Pi/2:+Pi/2] and lampda[-Pi:Pi]
	 * */
	public double lampda2Kavrayskiy(double lampdaRad, double phiRad){
		double lampda = (3*lampdaRad/(2*Math.PI))*Math.sqrt((Math.PI*Math.PI/3)-(phiRad*phiRad));
//		double test = (2*Math.PI*lampda) / (3*Math.sqrt((Math.PI*Math.PI/3)-(phiRad*phiRad))) ;
		return lampda;
	}
	public double lampdaFromKavrayskiy(double lampdaRad, double phiRad){
//		double lampda = (3*lampdaRad/(2*Math.PI))*Math.sqrt((Math.PI*Math.PI/3)-(phiRad*phiRad));
		double lampda = (2*Math.PI*lampdaRad) / (3*Math.sqrt((Math.PI*Math.PI/3)-(phiRad*phiRad))) ;
		return lampda;
	}
	
	/*
	 * Orthographic projection of angle values from certain viewpoint
	 * Expects values for phi[-Pi/2:+Pi/2] and lampda[-Pi:Pi]
	 * */
	private Pair<Double> orthoProjOfLP(double lampdaRad, double phiRad){
		Pair<Double> pos;
		double xPos=0, yPos=0;	
//		double lampda0 = this.origCoordinates.getFirst(); //Math.PI;
//		double phi0 = this.origCoordinates.getSecond(); //Math.PI/2;
		double lampda0 = this.centerOfProjection.getFirst(); //Math.PI;
		double phi0 = this.centerOfProjection.getSecond(); //Math.PI/2;
		lampda0 -= Math.PI;
		phi0 -= (Math.PI/2);		
		// values of lampdaRad and phiRad should lie between -Pi:+PI and -PI/2:+PI/2
		xPos = Math.cos(phiRad)*Math.sin(lampdaRad-lampda0);
		yPos = (Math.cos(phi0)*Math.sin(phiRad)) - (Math.sin(phi0)*Math.cos(phiRad)*Math.cos(lampdaRad-lampda0));
		pos = new Pair<Double>(xPos, yPos);
		return pos;
	}
	
	private Pair<Double> lpFromOrthoProj(double x, double y){
		Pair<Double> lp;
		double l=0, p=0;
//		boolean front = true;
//		if (x>=2*this.rSphere){
//			front = false;
//			x -= (2*this.rSphere);
//		}
//		x -= this.rSphere;
//		y -= this.rSphere;
		double dist = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2));
		double c = Math.asin(dist);
		double l0 = this.centerOfProjection.getFirst() - Math.PI; 
		double p0 = this.centerOfProjection.getSecond() - (Math.PI/2);
		p = Math.asin( Math.cos(c)*Math.sin(p0) + (y*Math.sin(c)*Math.cos(p0)/dist) );
		l = l0 + Math.atan( x*Math.sin(c) / ( (dist*Math.cos(p0)*Math.cos(c)) - (y*Math.sin(p0)*Math.sin(c)) ) );
		l += Math.PI;
		p += (Math.PI/2);
//		if (!front){
//			
//		}			
		lp = new Pair<Double>(l, p);
		return lp;
	}
	
	/*
	 * Transforms Radian values to Screen-position values
	 * Expects values for phi[0:Pi] and lampda[0:2Pi]
	 * */
	private Pair<Double> getScreenPosFromRad(double lampdaRad, double phiRad){
		Pair<Double> pos; // = new Pair(0,0);
		double xPos=0, yPos=0;	
		// 0:lampdaRad:2PI 0:phiRad:PI/2 
		
		if (this.mapProjType==kavrayskiyMapProj){
			phiRad = translateYCoordRespective2Orig(phiRad);
			lampdaRad = translateXCoordRespective2Orig(lampdaRad);
			lampdaRad -= Math.PI;
			phiRad -= (Math.PI/2);
			
			yPos = phiRad +(Math.PI/2);
			xPos = lampda2Kavrayskiy(lampdaRad, phiRad);		
			xPos = ( (Math.PI+xPos) *this.voxelsize )+this.border;
			yPos = (yPos*this.voxelsize) +this.border;				
		}
		else if (this.mapProjType==cylindricalMapProj){
			xPos = translateXPixelCoordRespective2Orig( lampdaRad *this.voxelsize )+this.border;
			yPos = translateYPixelCoordRespective2Orig(phiRad*this.voxelsize) +this.border;
		}
		else if (this.mapProjType==azimuthalMapProj){
//			double R = this.g2dSize.getHeight() / 2;
			// angles negative and positive values
			phiRad -= (Math.PI/2);
			lampdaRad -= Math.PI;	
			xPos = orthoProjOfLP(lampdaRad, phiRad).getFirst();
			yPos = orthoProjOfLP(lampdaRad, phiRad).getSecond();
			xPos = (xPos*this.voxelsize*voxelsizeFactor)+this.border;
			yPos = (yPos*this.voxelsize*voxelsizeFactor)+this.border;			
//			xPos = R * Math.cos(Math.abs(phiRad)) * Math.cos(lampdaRad);
//			yPos = R * Math.cos(Math.abs(phiRad)) * Math.sin(lampdaRad);
			xPos += this.rSphere; //R;
			yPos += this.rSphere; //R;
		}		
		pos = new Pair<Double>(xPos, yPos);
		return pos;
	}
	
	private boolean isOnFrontView(double lampda, double phi){
		double dist = distPointsOnSphere(lampda, phi, this.centerOfProjection.getFirst(), this.centerOfProjection.getSecond());
		if (dist>this.maxDistPoints)
			return false;
		else
			return true;
	}
	
	/**
	 * Returns the corresponding phi-lampda values in the sphoxelView given screen
	 * coordinates
	 */
	private Pair<Double> screen2spherCoord(Point point){
		
		Pair<Double> doubleP; // = new Pair<Double>((double)point.x/this.voxelsize , (double)point.y/this.voxelsize); 
		if (this.mapProjType==azimuthalMapProj){
			double x = point.x, y = point.y;
			if (x>(2*this.rSphere))
				x -= (2*this.rSphere);
			x -= this.rSphere;
			y -= this.rSphere; //R;
			doubleP = new Pair<Double>(x/(this.voxelsize*voxelsizeFactor) , y/(this.voxelsize*voxelsizeFactor)); 
		}
		else 
			doubleP = new Pair<Double>((double)point.x/this.voxelsize , (double)point.y/this.voxelsize); 
		
		return doubleP;
	}
	
	@SuppressWarnings("unused")
	private Pair<Double> translateCoordRespective2Orig(Pair<Double> pair){
		Pair<Double> transl;
		double xPos = translateXCoordRespective2Orig(pair.getFirst());
		double yPos = translateYCoordRespective2Orig(pair.getSecond());
		transl = new Pair<Double>(xPos, yPos);
		return transl;
	}
	
	private double translateXCoordRespective2Orig(double x){
		double dx = (this.origCoordinates.getFirst() - Math.PI);
		double xPos = x + dx;
		if (xPos > 2*Math.PI)
			xPos -= (2*Math.PI);
		else if (xPos < 0)
			xPos += (2*Math.PI);
		return xPos;
	}
	
	private double translateYCoordRespective2Orig(double y){
		double dy = (this.origCoordinates.getSecond() - (Math.PI/2));
		double yPos = y + dy;
//		yPos = (int)(yPos*1000)/1000;
//		if (yPos > Math.PI)
//			yPos -= Math.PI;
//		else if (yPos < 0)
//			yPos += Math.PI;
		return yPos;
	}
	
	private double translateXPixelCoordRespective2Orig(double x){
		double dx = (this.origCoordinates.getFirst() - Math.PI) * this.voxelsize;
		double xPos = x + dx;
		xPos = (int)(xPos*1000)/1000;
		if (xPos > 2*Math.PI* this.voxelsize)
			xPos -= (2*Math.PI* this.voxelsize);
		else if (xPos < 0)
			xPos += (2*Math.PI* this.voxelsize);
		return xPos;
	}
	
	private double translateYPixelCoordRespective2Orig(double y){
		double dy = (this.origCoordinates.getSecond() - (Math.PI/2)) * this.voxelsize;
		double yPos = y + dy;
//		yPos = (int)(yPos*1000)/1000;
//		if (yPos > Math.PI* this.voxelsize)
//			yPos -= Math.PI* this.voxelsize;
//		else if (yPos < 0)
//			yPos += Math.PI* this.voxelsize;
		return yPos;
	}
	
	/**
	 * Deleted selected range from list.
	 */
	public void deleteSelectedRange(){
		
		this.lampdaRanges.removeElementAt(this.chosenSelection);
		this.phiRanges.removeElementAt(this.chosenSelection);
		this.selContacts.removeElementAt(this.chosenSelection);	
	}
	
	private int checkForSelectedRanges(){
		int index = -1;	
		int count = 0;

		if (this.selContacts.size()>0){
			Iterator<Pair<Integer>> itrS = this.selContacts.iterator();
			while (itrS.hasNext()){
				Pair<Integer> sel = itrS.next();
				if (sel.getFirst()==this.tmpSelContact.getFirst() && sel.getSecond()==this.tmpSelContact.getSecond()){
					index = count; 
					break;
				}
				count++;
			}			
		}
	
		return index;
	}
	
	/**
	 * Update tmpContact with the contacts contained in the rectangle given by
	 * the upperLeft and lowerRight.
	 */
	private void squareSelect(){
		Pair<Double> upperLeft = screen2spherCoord(mousePressedPos);
		Pair<Double> lowerRight = screen2spherCoord(mouseDraggingPos);
		double xPosUL, yPosUL, xPosLR, yPosLR;
		// we reset the tmpContacts list so every new mouse selection starts
		// from a blank list
//		tmpContacts = new IntPairSet();
		
		xPosUL = upperLeft.getFirst();
		yPosUL = upperLeft.getSecond();
		xPosLR = lowerRight.getFirst();
		yPosLR = lowerRight.getSecond();
		if (this.mapProjType==kavrayskiyMapProj){
			xPosUL = lampdaFromKavrayskiy(xPosUL-Math.PI, yPosUL-(Math.PI/2)) + Math.PI;
			xPosLR = lampdaFromKavrayskiy(xPosLR-Math.PI, yPosLR-(Math.PI/2)) + Math.PI;
		}
		xPosUL -= (this.origCoordinates.getFirst()-Math.PI);
		xPosLR -= (this.origCoordinates.getFirst()-Math.PI);

		double pmin = Math.min(xPosUL, xPosLR);
		double tmin = Math.min(yPosUL, yPosLR);
		double pmax = Math.max(xPosUL, xPosLR);
		double tmax = Math.max(yPosUL, yPosLR);

//		System.out.println("squareSelect "+pmin+"-"+pmax+" , "+tmin+"-"+tmax);
//		boolean valid = true;
//		if (angleWithinValidRange(pmin, tmin) || angleWithinValidRange(pmax, tmax))
//			valid = false;
//		if (!valid){
//			pmin = 0;
//			tmin = 0;
//			pmax = 0;
//			tmax = 0;
//		}
//		if (valid)
		{
			this.tmpLampdaRange = new Pair<Double>(pmin, pmax);
			this.tmpPhiRange = new Pair<Double>(tmin, tmax);
//			this.tmpSelContact = new Pair<Integer>(this.AAStr.indexOf(this.iRes),this.AAStr.indexOf(this.jRes));
			this.tmpSelContact = new Pair<Integer>(this.iNum,this.jNum);
			
			this.tmpAngleComb = new Pair<Double>(pmax, tmax);
		}
	}
	
//	/** Called by ResidueRuler to enable display of ruler "crosshair" */	
//	public void showRulerCrosshair() {
//		showRulerCrosshair = true;
//	}
//	/** Called by ResidueRuler to switch off showing ruler coordinates */
//	public void hideRulerCoordinate() {
//		showRulerCoord = false;
//	}
	
	private void updateDimensions(){
		int xDimBU = this.xDim;
		int yDimBU = this.yDim;
		this.screenSize = this.contactView.getScreenSize();
		
		if ((float)(this.getSize().height)/(float)(this.getSize().width) > g2dRatio){
			this.xDim = this.getSize().width;
			this.yDim = (int) (this.xDim*g2dRatio);
		}
		else{
			this.yDim = this.getSize().height;
			this.xDim = (int) (this.yDim/g2dRatio);			
		}
		this.g2dSize.setSize(new Dimension(this.xDim, this.yDim));
		this.setSize(this.g2dSize);
		
		this.rSphere = g2dSize.getHeight()/2;
		this.maxDistPoints = distPointsOnSphere(0, 0, 0, Math.PI/2);
		
		// update pixel dimensions
		this.pixelWidth = (float)(this.xDim-2*this.border)/(float)(2*this.numSteps) ;
		this.pixelHeight =  this.pixelWidth; 
		this.voxelsize = (float) ((float)(this.numSteps*this.pixelWidth)/Math.PI);
		
		if (this.xDim!=xDimBU || this.yDim!=yDimBU)
			updateScreenBuffer();
		
	}

	/**
	 * Main method to draw the component on screen. This method is called each
	 * time the component has to be (re) drawn on screen. It is called
	 * automatically by Swing or by explicitly calling cmpane.repaint().
	 */
	@Override
	protected synchronized void paintComponent(Graphics g) {
//		System.out.println("paintComponent");
		
		Graphics2D g2d = (Graphics2D) g.create();		
//		updateDimensions();
		
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		// draw screen buffer
		if(screenBuffer==null) {
			screenBuffer = new ScreenBuffer(this);
			updateScreenBuffer();
		} else {
			
			g2d.drawImage(screenBuffer.getImage(),0,0,this);
		}
		// drawing selection rectangle if dragging mouse and showing temp
		// selection in red (tmpContacts)		
		drawOccupiedAngleRanges(g2d);
		drawSelectedAngleRange(g2d);
		drawCrosshair(g2d);
		updateAnglePanel();
	}
	
	/**
	 * Repaint the screen buffer because something in the underlying data has
	 * changed.
	 */
	private synchronized void updateScreenBuffer() {

		if(screenBuffer == null) {
			screenBuffer = new ScreenBuffer(this);
		}
		screenBuffer.clear();
		Graphics2D g2d = screenBuffer.getGraphics();

		// paint background
		int bgSizeX = Math.max(this.g2dSize.width, this.getWidth());		// fill background
															// even if window is
															// not square
		int bgSizeY = Math.max(this.g2dSize.height, this.getHeight());	// fill background
															// even of window is
															// not square
				
		g2d.setColor(backgroundColor);
		if (isOpaque()) {
			g2d.fillRect(0, 0, bgSizeX, bgSizeY);
		}
		
		// sphoxel background
		if (this.contactView.getGUIState().getShowSphoxelBG()){
//			drawSphoxels(g2d);
			drawSphoxelMap(g2d);
		}
		// neighbourhood-string-traces
		if (this.contactView.getGUIState().getShowNBHStraces()){
			drawNBHSTraces(g2d);
		}
		if (this.contactView.getGUIState().getShowLongitudes()){
			drawLongitudes(g2d);
		}
		if (this.contactView.getGUIState().getShowLatitudes()){
			drawLatitudes(g2d);
		}
		if (this.contactView.getGUIState().getShowLongLatCentre()){
			drawLongLatCentre(g2d);
		}				

		repaint();
	}
	
	/*---------------------------- setters and getters -----------------------------*/
	
	public String getNbhString() {
		return nbhString;
	}
	public String getOrigNBHString() {
		return origNBHString;
	}
	public void setNbhString(String nbhString) {
		this.nbhString = nbhString;
		this.nbhStringL = "%";
		this.nbhsRes = new char[this.nbhString.length()];
		int count = 0;
		this.nbhsRes = new char[this.nbhString.length()];
		for (int i=0; i<this.nbhString.length(); i++){
			this.nbhStringL += this.nbhString.charAt(i);
			this.nbhStringL += "%";
			this.nbhsRes[count] = this.nbhString.charAt(i);
			count++;
		}
		System.out.println(" Actual NBHS: "+this.nbhString+"-->"+this.nbhStringL);	
	}

	
	public void setStatusBar(ContactStatusBar statusBar) {
		this.contStatBar = statusBar;		
	}
	public ContactStatusBar getSatusBar(){
		return this.contStatBar;
	}
	
	public Dimension getScreenSize(){
		return this.screenSize;
	}
	
	public Dimension getPreferredSize(){
		return this.g2dSize;
	}
	
	public Dimension getPanelSize(){
		return this.g2dSize;
	}

	public void setRadiusRangesFixed(boolean radiusRangesFixed) {
		this.radiusRangesFixed = radiusRangesFixed;
	}

	public boolean isRadiusRangesFixed() {
		return radiusRangesFixed;
	}
	
	public void setMinR(float r){
		this.minr = r;
	}
	public void setMaxR(float r){
		this.maxr = r;
	}
	public void setNumSteps(int num){
		this.numSteps = num;
		this.resol = 180/(float)this.numSteps;
		this.deltaRad = Math.PI/this.numSteps;
	}
	public void setResol(float res){
		this.numSteps = (int) (180/res);
		this.resol = 180/(float)this.numSteps;
		System.out.println("numSteps="+this.numSteps+"  --> resol="+this.resol);
	}	
	public Pair<Double> getOrigCoord(){
		return this.origCoordinates;
	}
	public float getVoxelsize(){
		return this.voxelsize;
	}

	public boolean isDiffSStype() {
		return diffSStype;
	}
	public void setDiffSStype(boolean diffSStype) {
		this.diffSStype = diffSStype;
		if (!this.diffSStype)
			this.iSSType = CMPdb_sphoxel.AnySStype;
	}

	public boolean isRemoveOutliers() {
		return removeOutliers;
	}

	public void setRemoveOutliers(boolean removeOutliers) {
		this.removeOutliers = removeOutliers;
	}

	public double getMinAllowedRat() {
		return minAllowedRat;
	}

	public void setMinAllowedRat(double minAllowedRat) {
		this.minAllowedRat = minAllowedRat;
		updateScreenBuffer();
	}

	public double getMaxAllowedRat() {
		return maxAllowedRat;
	}

	public void setMaxAllowedRat(double maxAllowedRat) {
		this.maxAllowedRat = maxAllowedRat;
		updateScreenBuffer();
	}


	public int getChosenColourScale() {
		return chosenColourScale;
	}

	public void setChosenColourScale(int chosenColourScale) {
		this.chosenColourScale = chosenColourScale;
		// show actual colour scale
		updateScreenBuffer();
	}

	public int getMapProjType() {
		return mapProjType;
	}

	public void setMapProjType(int mapProjType) {
		this.mapProjType = mapProjType;
		if (this.mapProjType!=azimuthalMapProj)
			this.origCoordinates = new Pair<Double>(this.origCoordinates.getFirst(), Math.PI/2);
		else{
			this.origCoordinates = new Pair<Double>(this.origCoordinates.getFirst(), this.centerOfProjection.getSecond());
//			this.origCoordinates = this.centerOfProjection;
			this.centerOfProjection = this.origCoordinates;
		}
		this.updateScreenBuffer();
	}

	public void setDeltaOffSetX(double deltaOffSetX) {
		this.deltaOffSetX = deltaOffSetX;
	}

	public double getDeltaOffSetX() {
		updateDimensions();
		this.deltaOffSetX = getOffSet(0.0, 0.0);
		return deltaOffSetX;
	}
	
	public void setDeltaOffSetXEnd(double deltaOffSetXEnd) {
		this.deltaOffSetXEnd = deltaOffSetXEnd;
	}

	public double getDeltaOffSetXEnd() {
		updateDimensions();
		this.deltaOffSetXEnd = getOffSet(2*Math.PI, 0.0);
		return deltaOffSetXEnd;
	}
	
	public void setDeltaOffSetXCenter(double deltaOffSetXC) {
		this.deltaOffSetXCenter = deltaOffSetXC;
	}

	public double getDeltaOffSetXCenter() {
		updateDimensions();
		// this.origCoordinates = new Pair<Float>((float)(Math.PI), (float)(Math.PI/2));
		this.deltaOffSetXCenter = getOffSet(this.origCoordinates.getFirst(), 0.0);
		return deltaOffSetXCenter;
	}

	public double getMinRatio() {
		return minRatio;
	}

	public void setMinRatio(double minRatio) {
		this.minRatio = minRatio;
	}

	public double getMaxRatio() {
		return maxRatio;
	}

	public void setMaxRatio(double maxRatio) {
		this.maxRatio = maxRatio;
	}	

	public int[] getHistWholeAngleRange() {
		return histWholeAngleRange;
	}
	public Vector<int[]> getHist4SelAngleRange() {
		return hist4SelAngleRange;
	}
	public void setShowHist(boolean showHist) {
		this.showHist = showHist;
	}
	public boolean isShowHist() {
		return showHist;
	}	
	public void setHistType(int type){
		this.histType = type;
	}
	public int getHistType(){
		return this.histType;
	}
	public Vector<double[]> getMinMaxAverT4SelAngleRange() {
		return minMaxAverT4SelAngleRange;
	}
	public void setMinMaxAverT4SelAngleRange(
			Vector<double[]> minMaxAverT4SelAngleRange) {
		this.minMaxAverT4SelAngleRange = minMaxAverT4SelAngleRange;
	}
	public int[][] getHistTWholeAngleRange() {
		return histTWholeAngleRange;
	}	
	public void setHistTWholeAngleRange(int[][] histTWholeAngleRange) {
		this.histTWholeAngleRange = histTWholeAngleRange;
	}
	public Vector<int[][]> getHistT4SelAngleRange() {
		return histT4SelAngleRange;
	}
	public void setHistT4SelAngleRange(Vector<int[][]> histT4SelAngleRange) {
		this.histT4SelAngleRange = histT4SelAngleRange;
	}

	/* ChosenSelection has value [0:numSelections-1].
	 * If right click outside of any selection, chosen selection = numSelections.
	 *  */
	public int getChosenSelection() {
		return chosenSelection;
	}
	public void setChosenSelection(int chosenSelection) {
		this.chosenSelection = chosenSelection;
	}
	
	public String[] getChosenPTRange(){
//		double[] range;
		String[] rangeS;
		Iterator<Pair<Double>> itrP = this.lampdaRanges.iterator();
		Iterator<Pair<Double>> itrT = this.phiRanges.iterator();
		int count = 0;
		Pair<Double> lampda=null, phi=null;
		while (count<=this.chosenSelection && itrP.hasNext()){
			lampda = itrP.next();
			phi = itrT.next();
			count++;
		}
		double p1 = (double)(Math.round(lampda.getFirst()*100))/100;
		double p2 = (double)(Math.round(lampda.getSecond()*100))/100;
		double t1 = (double)(Math.round(phi.getFirst()*100))/100;
		double t2 = (double)(Math.round(phi.getSecond()*100))/100; //(float)(Math.round(this.tmpAngleComb.getFirst()*100))/100;
		rangeS = new String[]{String.valueOf(p1), String.valueOf(p2), String.valueOf(t1), String.valueOf(t2)};
//		range = new double[]{lampda.getFirst(), lampda.getSecond(), phi.getFirst(), phi.getSecond()};
		return rangeS;
	}

	
//	/**
//	 * Sets the output size and updates the ratio and contact square size. This
//	 * will affect all drawing operations. Used by print() method to change the
//	 * output size to the size of the paper and back.
//	 */
//	protected void setOutputSize(int size) {
//		outputSize = size;
//		ratio = (double) outputSize/numSteps;		// scale factor, = size
////		// of one pixel
////		this.pixelWidth =  (float) ratio; 			// the size of the
////		this.pixelHeight = (float) ratio;
////		// square representing a sphoxel
//	}
	
	private void showPopup(MouseEvent e) {
//		this.rightClickAngle = screen2cm(new Point(e.getX(), e.getY()));
		// determine chosen selection
		Point point = e.getPoint();
		Pair<Double> pos = screen2spherCoord(point);
		double xPos = pos.getFirst();
		double yPos = pos.getSecond();
		
		if (this.mapProjType==kavrayskiyMapProj){
			xPos = lampdaFromKavrayskiy(xPos-Math.PI, yPos-(Math.PI/2)) + Math.PI;
		}
		xPos -= (this.origCoordinates.getFirst()-Math.PI);		
		this.rightClickAngle = new Pair<Double>(xPos, yPos);
		
		Iterator<Pair<Double>> itrP = this.lampdaRanges.iterator();
		Iterator<Pair<Double>> itrT = this.phiRanges.iterator();
		int count = 0;
		this.chosenSelection = this.lampdaRanges.size();
		while (itrP.hasNext()){
			Pair<Double> lampda = itrP.next();
			Pair<Double> phi = itrT.next();
			if (rightClickAngle.getFirst()>=lampda.getFirst() && rightClickAngle.getFirst()<=lampda.getSecond() 
					&& rightClickAngle.getSecond()>=phi.getFirst() && rightClickAngle.getSecond()<=phi.getSecond()){
				this.chosenSelection = count;
				break;
			}
			count++;
		}
		System.out.println("chosen selection= "+this.chosenSelection);
//		if (this.lampdaRanges.size()>0){			
//		}
//		else {			
//		}		
		if (this.chosenSelection != this.lampdaRanges.size())
			this.contactView.popup.show(e.getComponent(), e.getX(), e.getY());
	}
	
	/*---------------------------- mouse events -----------------------------*/

	public void mousePressed(MouseEvent evt) {
		// This is called when the user presses the mouse anywhere
		// in the frame
		
		System.out.println("mousePressed");

		lastMouseButtonPressed = evt.getButton();
		mousePressedPos = evt.getPoint();
		if(lastMouseButtonPressed == MouseEvent.BUTTON2) dragging = false;
		
		if (evt.isPopupTrigger()) {
			showPopup(evt);
			return;
		}
	}

	public void mouseReleased(MouseEvent evt) {

		// Called whenever the user releases the mouse button.
		// TODO: Move much of this to MouseClicked and pull up Contact cont = screen2spherCoord...
		if (evt.isPopupTrigger()) {
			showPopup(evt);
			dragging = false;
			return;
		}
		// only if release after left click (BUTTON1)
		if (evt.getButton()==MouseEvent.BUTTON1) {

			switch (contactView.getGUIState().getSelectionMode()) {
			case RECT:
				if (this.mapProjType!=azimuthalMapProj){					
					if (dragging){
	//					squareSelect();
						boolean valid = true; // = angleWithinValidRange(this.tmpSelContact.getFirst(), this.tmpSelContact.getSecond());
						if (!angleWithinValidRange(tmpLampdaRange.getFirst(), tmpPhiRange.getFirst()) 
								|| !angleWithinValidRange(tmpLampdaRange.getSecond(), tmpPhiRange.getSecond()))
							valid = false;
						if (valid){
							int index = checkForSelectedRanges();
							if (index>-1){
								this.lampdaRanges.removeElementAt(index);
								this.phiRanges.removeElementAt(index);
								this.selContacts.removeElementAt(index);												
							}
							this.lampdaRanges.add(this.tmpLampdaRange);
							this.phiRanges.add(this.tmpPhiRange);
							this.selContacts.add(this.tmpSelContact);
							updateAngleRange();
							System.out.println("mouseReleased "+this.tmpLampdaRange.getFirst()+"-"+this.tmpLampdaRange.getSecond()
									+" , "+this.tmpPhiRange.getFirst()+"-"+this.tmpPhiRange.getSecond());						
						}
					}				
					dragging = false;
	//				this.repaint();				
				}
				
				return;
				
			case CLUSTER:
//				// resets selContacts when clicking mouse

//				this.repaint();
				return;
				
			case PAN:
				if (this.mapProjType != azimuthalMapProj){
					Pair<Double> pos = screen2spherCoord(mousePos);
					double xPos = pos.getFirst();
					double yPos = Math.PI/2; //pos.getSecond();
					if (this.mapProjType==kavrayskiyMapProj){
						xPos = lampdaFromKavrayskiy(xPos-Math.PI, yPos-(Math.PI/2)) + Math.PI;
					}
					this.origCoordinates = new Pair<Double>(xPos, yPos);
//					this.centerOfProjection = new Pair<Double>(this.origCoordinates.getFirst()-Math.PI, this.origCoordinates.getSecond()-Math.PI/2);
//					this.centerOfProjection = this.origCoordinates;
					
					this.contactView.lampdaRuler.repaint();
					this.contactView.phiRuler.repaint();
					updateScreenBuffer();
					System.out.println("tmpAngleComb: "+this.tmpAngleComb.getFirst()+","+this.tmpAngleComb.getSecond());
					System.out.println("OrigCoord: "+this.origCoordinates.getFirst()+","+this.origCoordinates.getSecond());
					System.out.println("CentreCoord: "+this.centerOfProjection.getFirst()+","+this.centerOfProjection.getSecond());
//					this.repaint();					
				}
				return;
			}
			
			this.tmpLampdaRange = new Pair<Double>(0.0, 0.0);
			this.tmpPhiRange = new Pair<Double>(0.0, 0.0);
		}
	}

	public void mouseDragged(MouseEvent evt) {
//		System.out.println("mouseDragged");

		// Called whenever the user moves the mouse
		// while a mouse button is held down.

		mouseMoved(evt); // TODO is this necessary? I tried getting rid of it
		// but wasn't quite working

		if(lastMouseButtonPressed == MouseEvent.BUTTON1) {
			dragging = true;
			mouseDraggingPos = evt.getPoint();
			switch (contactView.getGUIState().getSelectionMode()) {
			case RECT:
				if (this.mapProjType!=azimuthalMapProj)
					squareSelect();				
				break;
			case PAN:
				if (this.mapProjType != azimuthalMapProj){
					Pair<Double> pos = screen2spherCoord(mousePos);
					double xPos = pos.getFirst();
					double yPos = Math.PI/2; //pos.getSecond();
					double fac = 2.5*2*Math.PI/360; //(Math.PI / 20);
					if (this.mapProjType==kavrayskiyMapProj){
						xPos = lampdaFromKavrayskiy(xPos-Math.PI, yPos-(Math.PI/2)) + Math.PI;
					}
					xPos = xPos/fac;
					xPos = Math.round(xPos);
					xPos = xPos*fac;
//					xPos = Math.round/(xPos/fac) * fac;
					this.origCoordinates = new Pair<Double>(xPos, yPos);
//					this.centerOfProjection = new Pair<Double>(Math.PI-this.origCoordinates.getFirst(), Math.PI/2-this.origCoordinates.getSecond());
//					this.centerOfProjection = new Pair<Double>(this.origCoordinates.getFirst()-Math.PI, this.origCoordinates.getSecond()-Math.PI/2);
//					this.centerOfProjection = this.origCoordinates;
					this.contactView.lampdaRuler.repaint();
					this.contactView.phiRuler.repaint();
					updateScreenBuffer();					
				}
				break;
			}
		}
	} 

	public void mouseEntered(MouseEvent evt) { 
		mouseIn = true;
	}

	public void mouseExited(MouseEvent evt) {
//		System.out.println("mouseExited");

		mouseIn = false;
//		this.repaint();
	}

	public void mouseClicked(MouseEvent evt) {
	}
	
	public void mouseMoved(MouseEvent evt) {
//		System.out.println("mouseMoved");
		double xPos, yPos;
		double l, p;
		mousePos = evt.getPoint();
		
		Pair<Double> pos = screen2spherCoord(mousePos);
		xPos = pos.getFirst();
		yPos = pos.getSecond();
		if (this.mapProjType==kavrayskiyMapProj){
			l = lampdaFromKavrayskiy(xPos-Math.PI, yPos-(Math.PI/2)) + Math.PI;
			p = yPos;
		}
		else if (this.mapProjType==azimuthalMapProj){
			l = lpFromOrthoProj(xPos, yPos).getFirst();
			p = lpFromOrthoProj(xPos, yPos).getSecond();
		}
		else{
			l = xPos;
			p = yPos;
		}
		boolean valid = angleWithinValidRange(l, p);
		l -= (this.origCoordinates.getFirst()-Math.PI);
		if (l<0)
			l += (2*Math.PI);
		else if (l>2*Math.PI)
			l -= (2*Math.PI);
		
//		this.tmpAngleComb = new Pair<Double>(pos.getFirst(), pos.getSecond());
//		boolean valid = angleWithinValidRange(xPos, yPos);
//		if (this.mapProjType==azimuthalMapProj)
//			valid = false;
		
		if (valid)
			this.tmpAngleComb = new Pair<Double>(l, p);
		else 
			this.tmpAngleComb = new Pair<Double>(0.0, 0.0);
//		this.tmpAngleComb = new Pair<Double>((float)xPos, (float)yPos);
		this.repaint();
	}

	public void keyPressed(KeyEvent e) {
		// TODO Auto-generated method stub
//		System.out.println("keyPressed");
		
	}

	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub
//		System.out.println("keyReleased");
		//-- Process arrow "virtual" keys
    	double first = this.origCoordinates.getFirst();
    	double second = this.origCoordinates.getSecond();
    	double fac = 2*(2.5*2*Math.PI/360); //(Math.PI / 20);
        switch (e.getKeyCode()) {
            case KeyEvent.VK_LEFT : first-=fac; break;
            case KeyEvent.VK_RIGHT: first+=fac; break;
            case KeyEvent.VK_UP   : second-=fac;   break;
            case KeyEvent.VK_DOWN : second+=fac; break;
        }
        if (first<0)
        	first += 2*Math.PI;
        if (first>2*Math.PI)
        	first -= 2*Math.PI;
        if (second<0)
        	second += Math.PI;
        if (second>Math.PI)
        	second -= Math.PI;
        
        if (this.mapProjType!=azimuthalMapProj)
        	second = (float) (Math.PI/2);   // --> panning just horizontally     	
		this.origCoordinates = new Pair<Double>(first,second); // this.tmpAngleComb;
		if (this.mapProjType==azimuthalMapProj){
//			this.centerOfProjection = this.origCoordinates;
			this.centerOfProjection = new Pair<Double>(2*Math.PI-this.origCoordinates.getFirst(), Math.PI-this.origCoordinates.getSecond());
		}

		System.out.println("OrigCoord: "+this.origCoordinates.getFirst()+","+this.origCoordinates.getSecond());
		System.out.println("CentreCoord: "+this.centerOfProjection.getFirst()+","+this.centerOfProjection.getSecond());
		
		this.contactView.lampdaRuler.repaint();
		this.contactView.phiRuler.repaint();	
		updateScreenBuffer();
	}

	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub	
//		System.out.println("keyTyped");	
	}

	@Override
	public void componentHidden(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void componentMoved(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void componentResized(ComponentEvent e) {
		// TODO Auto-generated method stub
		updateDimensions();
	}

	@Override
	public void componentShown(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}
}
