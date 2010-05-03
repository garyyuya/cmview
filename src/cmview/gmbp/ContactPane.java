package cmview.gmbp;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JPanel;

import owl.core.structure.features.SecStrucElement;
import owl.core.structure.graphs.RIGNbhood;
import owl.core.structure.graphs.RIGNode;

import owl.gmbp.CMPdb_nbhString_traces;
import owl.gmbp.CMPdb_sphoxel;
import owl.gmbp.CSVhandler;

import cmview.ContactMapPane;
import cmview.Start;
import cmview.datasources.Model;
import edu.uci.ics.jung.graph.util.Pair;

public class ContactPane extends JPanel implements MouseListener, MouseMotionListener, KeyListener{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	protected static final Dimension defaultDim = new Dimension(1600, 800);
	protected static final float g2dRatio = 0.5f; // H/W
	protected static final double defaultMinAllowedRat = -3;
	protected static final double defaultMaxAllowedRat = 1;
		
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
	
	// drawing colors (being set in the constructor)
	private Color backgroundColor;	  		// background color
	private Color squareSelColor;	  		// color of selection rectangle
	private Color crosshairColor;     		// color of crosshair	
	private Color selAngleRangeColor;       // color for selected rectangles
	private Color longitudeColor;			// color for longitudes
	private Color latitudeColor;			// color for latitudes
	
	// selections 
//	private FloatPairSet phiRanges;			// permanent list of currently selected phi ranges
//	private FloatPairSet thetaRanges;       // permanent list of currently selected theta ranges
//	private IntPairSet selContacts;         // permanent list of currently selected and referred contacts
	private Vector<Pair<Float>> phiRanges;			// permanent list of currently selected phi ranges
	private Vector<Pair<Float>> thetaRanges;       // permanent list of currently selected theta ranges
	private Vector<Pair<Integer>> selContacts;         // permanent list of currently selected and referred contacts
	private Pair<Float> tmpPhiRange;
	private Pair<Float> tmpThetaRange;
	private Pair<Integer> tmpSelContact;
	private Pair<Float> tmpAngleComb;
	private Pair<Float> origCoordinates;    // actual coordinates of origin in angle range x:[0:2*PI] y:[0:PI]
	
	private RIGNode nodeI, nodeJ;
	
	// query data
	private char iRes='A', jRes='A';
	private char iSSType='H', jSSType='H';
	private boolean diffSStype=false;
	//	private String iResType="Ala", jResType="Ala";
	private int iNum=0, jNum=0;
	private String nbhString, nbhStringL;
	private String jAtom = "CA";
	private char[] nbhsRes;
	
	private String db = "bagler_all5p0";
	
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
	// NBHStraces-Data
	private CMPdb_nbhString_traces nbhsTraces;
	private Vector<float[]> nbhsNodes;
	private int[] numNodesPerLine;
	private int[] maxDistsJRes, minDistsJRes;
	private int numLines;
	
	// variables for representation (drawing methods)
	private int numSteps = CMPdb_sphoxel.defaultNumSteps; //CMPdb_sphoxel.defaultNumSteps;  // change later via interface
	private float resol = CMPdb_sphoxel.defaultResol; //CMPdb_sphoxel.getDefaultResol();
	private final int border = 0; //15;
	private final int yBorderThres = 0; //45;
	private float pixelWidth = 5*36/this.numSteps; // change in dependence of numSteps
	private float pixelHeight = this.pixelWidth; //5*36/this.numSteps;	
	private float voxelsize = (float) (this.numSteps*this.pixelWidth/Math.PI); //pixelWidth;

	private boolean removeOutliers = true;
	private double minAllowedRat = defaultMinAllowedRat;
	private double maxAllowedRat = defaultMaxAllowedRat;
	private int chosenColourScale = ContactStatusBar.BLUERED;
	
	private boolean paintCentralResidue = true;
	private boolean mapProjection = false;
	
	private int xDim=0;
	private int yDim=0;
	
	private final char[] aas = new char[]{'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}; // possible Residues
//	private final String AAStr = new String(aas); 
	
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
		this.setFocusable(true);
		this.addKeyListener(this);

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

		this.mousePos = new Point();
		this.mousePressedPos = new Point();
		this.mouseDraggingPos = new Point();
		
		// set default colors
		this.backgroundColor = Color.white;
		this.squareSelColor = Color.gray;	
		this.crosshairColor = Color.yellow;
		this.selAngleRangeColor = Color.black;
		this.longitudeColor = Color.black;
		this.latitudeColor = this.longitudeColor;
		
		this.dragging = false;
		this.selContacts = new Vector<Pair<Integer>>();
		this.phiRanges = new Vector<Pair<Float>>();
		this.thetaRanges = new Vector<Pair<Float>>();
		
		this.tmpPhiRange = new Pair<Float>(0.0f, 0.0f);
		this.tmpThetaRange = new Pair<Float>(0.0f, 0.0f);
		this.tmpAngleComb = new Pair<Float>(0.0f, 0.0f);
		this.origCoordinates = new Pair<Float>((float)(Math.PI), (float)(Math.PI/2));
		
//		setOutputSize(screenSize.width, screenSize.height);
//		System.out.println("outputsize= "+this.outputSize);
	}

	private void updateAngleRange(){
//		RIGEdge edge = this.mod.getGraph().getEdgeFromSerials(this.tmpSelContact.getFirst(), this.tmpSelContact.getSecond());		
//		edge.setPhiPsi(this.tmpPhiRange.getFirst(), this.tmpPhiRange.getSecond(), this.tmpThetaRange.getFirst(), this.tmpThetaRange.getSecond());		
		this.mod.getGmbp().setPhiRanges(this.phiRanges);
		this.mod.getGmbp().setThetaRanges(this.thetaRanges);
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
		this.jAtom = this.mod.edgeType;
		this.nbhString = nbhood.getNbString();
		this.nbhStringL = "%";
		int count = 0;
		this.nbhsRes = new char[this.nbhString.length()];
		for (int i=0; i<this.nbhString.length(); i++){
			this.nbhStringL += this.nbhString.charAt(i);
			this.nbhStringL += "%";
			this.nbhsRes[count] = this.nbhString.charAt(i);
			count++;
		}
		System.out.println(this.nbhString+"-->"+this.nbhStringL);	

	}
	
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
			
			String fn = Start.SPHOXEL_DIR; // "/Users/vehlow/Documents/workspace/outputFiles/sphoxelBG/";
			fn = fn+"sphoxelBG_"+this.iRes+"-"+String.valueOf(this.iSSType).toLowerCase()+"_"+this.jRes+"-"
				+String.valueOf(CMPdb_sphoxel.AnySStype).toLowerCase()+"_"+this.radiusPrefix+".csv";
			System.out.println("Filename= "+fn);
			fn = "/Users/vehlow/Documents/workspace/outputFiles/LogOddsScoresBayes_fromDB-bagler_all13p0_alledges_A-A_SStype-H_radius9.2-12.8_resol90.csv";
			
			CSVhandler csv = new CSVhandler();
			try {
				this.bayesRatios = csv.readCSVfile3Ddouble(fn);
				setNumSteps(this.bayesRatios.length);
				this.ratios = new double[this.bayesRatios.length][this.bayesRatios[0].length];
				for (int i=0; i<this.bayesRatios.length; i++){
					for (int j=0; j<this.bayesRatios[i].length; j++){
						this.ratios[i][j] = this.bayesRatios[i][j][0];
					}
				}
				for(int i=0;i<this.ratios.length;i++){ // dim for theta
					for(int j=0;j<this.ratios[i].length;j++){ // dim for phi					                 
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
		}
	}
	
	private void setTracesParam(){
		nbhsTraces.setDiffSSType(this.diffSStype);
		nbhsTraces.setSSType(this.iSSType);
	}
	
	private void calcNbhsTraces() throws SQLException{
		nbhsTraces.run();
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
	
	public void recalcTraces() throws SQLException{
		if (this.nbhsTraces.getDiffSSType()!=this.diffSStype || this.nbhsTraces.getSSType()!=this.iSSType){
			setTracesParam();
			calcNbhsTraces();
		}
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
	
	public GeneralPath nGon(Pair<Double> p1, Pair<Double> p2, Pair<Double> p3, Pair<Double> p4){
		GeneralPath shape = new GeneralPath();
		shape.moveTo(p1.getFirst(), p1.getSecond());
		shape.lineTo(p2.getFirst(), p2.getSecond());
		shape.lineTo(p3.getFirst(), p3.getSecond());
		shape.lineTo(p4.getFirst(), p4.getSecond());
		shape.lineTo(p1.getFirst(), p1.getSecond());
		return shape;
	}
	
	public GeneralPath nGon(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
		GeneralPath shape = new GeneralPath();
		shape.moveTo(x1, y1);
		shape.lineTo(x2, y2);
		shape.lineTo(x3, y3);
		shape.lineTo(x4, y4);
		shape.lineTo(x1, y1);
		return shape;
	}
	
	public double phiKavrayskiy(double phiRad, double thetaRad){
		double phi = (3*phiRad/(2*Math.PI))*Math.sqrt((Math.PI*Math.PI/3)-(thetaRad*thetaRad));
		return phi;
	}
	
	/*
	 * Transforms Radian values to Screen-position values
	 * Expects values for theta[0:Pi] and phi[0:2Pi]
	 * */
	public Pair<Double> getScreenPosFromRad(double phiRad, double thetaRad){
		Pair<Double> pos; // = new Pair(0,0);
		double xPos, yPos;	
		
		if (this.mapProjection){
			thetaRad = translateYCoordRespective2Orig(thetaRad);
			phiRad = translateXCoordRespective2Orig(phiRad);
			phiRad -= Math.PI;
			thetaRad -= (Math.PI/2);
			
			yPos = thetaRad +(Math.PI/2);
			xPos = phiKavrayskiy(phiRad, thetaRad);		
			xPos = ( (Math.PI+xPos) *this.voxelsize )+this.border;
			yPos = (yPos*this.voxelsize) +this.border;				
		}
		else {
			phiRad -= Math.PI;
			xPos = translateXPixelCoordRespective2Orig( (Math.PI+phiRad) *this.voxelsize )+this.border;
			yPos = translateYPixelCoordRespective2Orig(thetaRad*this.voxelsize) +this.border;
		}
		
		pos = new Pair<Double>(xPos, yPos);
		return pos;
	}
	
	public void drawSphoxelMap(Graphics2D g2d){
		Shape shape = null;
		double xPos1, yPos1, xPos2, yPos2, xPos3, yPos3, xPos4, yPos4;
		double val;
//		double thetaDeg, phiDeg; 
		double thetaRad, phiRad;
		Color col = Color.white; // = new Color(24,116,205,255);
		ColorScale scale = new ColorScale();
		double deltaRad = Math.PI/this.numSteps;
//		double deltaDeg = 180/this.numSteps;
		
		g2d.setBackground(backgroundColor);
		shape = new Rectangle2D.Float(0, 0, this.pixelWidth*this.ratios[0].length, this.pixelHeight*this.ratios.length);
		g2d.setColor(backgroundColor);
		g2d.draw(shape);
		g2d.fill(shape);
				
		// ----- color representation -----
		for(int i=0;i<ratios.length;i++){
//			if (i==ratios.length - 2){
//				int test = 0;
//				int ff = test;
//			}
			for(int j=0;j<ratios[i].length;j++){
//				thetaDeg = i*deltaDeg;
				thetaRad = (i*deltaRad);
//				phiDeg = 180 - (j*deltaDeg);
				phiRad = (j*deltaRad);
//				phiD = (j*deltaDeg);
//				phiRad = (j*deltaRad);
				
//				xPos1 = getScreenPosFromRad(phiRad, thetaRad).getFirst();
//				yPos1 = getScreenPosFromRad(phiRad, thetaRad).getSecond();
//				xPos2 = getScreenPosFromRad(phiRad+deltaRad, thetaRad).getFirst();
//				yPos2 = getScreenPosFromRad(phiRad+deltaRad, thetaRad).getSecond();	
//				thetaRad += deltaRad;
//				xPos3 = getScreenPosFromRad(phiRad+deltaRad, thetaRad).getFirst();
//				yPos3 = getScreenPosFromRad(phiRad+deltaRad, thetaRad).getSecond();
//				xPos4 = getScreenPosFromRad(phiRad, thetaRad).getFirst();
//				yPos4 = getScreenPosFromRad(phiRad, thetaRad).getSecond();
//				thetaRad -= deltaRad;
//				shape = nGon(xPos1, yPos1, xPos2, yPos2, xPos3, yPos3, xPos4, yPos4);
				
				if (this.mapProjection){
					thetaRad = translateYCoordRespective2Orig(thetaRad);
					phiRad = translateXCoordRespective2Orig(phiRad);
					phiRad -= Math.PI;
					thetaRad -= (Math.PI/2);
					
					yPos1 = thetaRad +(Math.PI/2);
					xPos1 = phiKavrayskiy(phiRad, thetaRad);
//					xPos1 = (3*phiRad/(2*Math.PI))*Math.sqrt((Math.PI*Math.PI/3)-(thetaRad*thetaRad));
					yPos2 = yPos1;
					xPos2 = phiKavrayskiy(phiRad+deltaRad, thetaRad);
//					xPos2 = (3*(phiRad+deltaRad)/(2*Math.PI))*Math.sqrt((Math.PI*Math.PI/3)-(thetaRad*thetaRad));	
					thetaRad += deltaRad;
					yPos3 = thetaRad +(Math.PI/2);
					xPos3 = phiKavrayskiy(phiRad+deltaRad, thetaRad);
//					xPos3 = (3*(phiRad+deltaRad)/(2*Math.PI))*Math.sqrt((Math.PI*Math.PI/3)-(thetaRad*thetaRad));	
					yPos4 = yPos3;
					xPos4 = phiKavrayskiy(phiRad, thetaRad);
//					xPos4 = (3*phiRad/(2*Math.PI))*Math.sqrt((Math.PI*Math.PI/3)-(thetaRad*thetaRad));
					
					xPos1 = ( (Math.PI+xPos1) *this.voxelsize )+this.border;
					yPos1 = (yPos1*this.voxelsize) +this.border;
					xPos2 = ( (Math.PI+xPos2) *this.voxelsize )+this.border;
					yPos2 = (yPos2*this.voxelsize) +this.border;
					xPos3 = ( (Math.PI+xPos3) *this.voxelsize )+this.border;
					yPos3 = (yPos3*this.voxelsize) +this.border;
					xPos4 = ( (Math.PI+xPos4) *this.voxelsize )+this.border;
					yPos4 = (yPos4*this.voxelsize) +this.border;

					shape = nGon(xPos1, yPos1, xPos2, yPos2, xPos3, yPos3, xPos4, yPos4);								
				}
				else {
//					xPos = translateXPixelCoordRespective2Orig( (Math.PI+phiRad) *this.voxelsize )+this.border;
//					yPos = translateYPixelCoordRespective2Orig(thetaRad*this.voxelsize) +this.border;				
					xPos1 = j*pixelWidth;
					yPos1 = i*pixelWidth;
					xPos1 = translateXPixelCoordRespective2Orig(xPos1) +this.border;
					yPos1 = translateYPixelCoordRespective2Orig(yPos1) +this.border;			
					shape = nGon(xPos1, yPos1, xPos1+pixelWidth, yPos1, xPos1+pixelWidth, yPos1+pixelHeight, xPos1, yPos1+pixelHeight);
//					shape = new Rectangle2D.Double(xPos1, yPos1, pixelWidth, pixelHeight);	
				}
				
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
				g2d.setColor(col);
				g2d.fill(shape);
//				if (xPos+pixelWidth > this.xDim)
//					xPos -= xDim;
//				if (yPos+pixelHeight > this.yDim)
//					yPos -= yDim;
//				shape = new Rectangle2D.Double(xPos, yPos, pixelWidth, pixelHeight);
//				g2d.draw(shape);
//				g2d.fill(shape);
			}
		}
	}
	
	public void drawSphoxels(Graphics2D g2d){
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
//				shape = nGon(xPos, yPos, xPos+pixelWidth, yPos, xPos+pixelWidth, yPos+pixelHeight, xPos, yPos+pixelHeight);
				
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
			col = (Color.magenta); break;
		case 1: // 'S'
			col = (Color.yellow); break;
		case 2: // 'O'
			col = (Color.cyan); break;
		}
		return col;
	}
	
	public void drawNBHSNode(Graphics2D g2d, boolean specialRes, float[] node){
		GeneralPath rhombus = null;
		Shape circle = null;
		int iNum, jNum, jResID, jSSType;
		double xPos, yPos;
		double thetaRad, phiRad;
		Font f = new Font("Dialog", Font.PLAIN, 12);
		float radius = 3.f;
		String nodeName;
		
		iNum = (int) node[1];
		jNum = (int) node[2];
		thetaRad = node[3];
		phiRad = node[4];
		jResID = (int) node[5];
		jSSType = (int) node[6];
		
		phiRad += Math.PI;
		
		xPos = getScreenPosFromRad(phiRad, thetaRad).getFirst();
		yPos = getScreenPosFromRad(phiRad, thetaRad).getSecond();
		
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
		nodeName = Character.toString(this.aas[jResID]) +" "+ String.valueOf(jNum-iNum);
		if (specialRes)
			g2d.setColor(Color.black);
		else
			g2d.setColor(new Color(70,70,70));
		g2d.drawString(nodeName, (float)(xPos+radius), (float)(yPos+radius+this.yBorderThres));
	}
	
	public void drawNBHSEdge(Graphics2D g2d, float[] node, int nodeID, int lineID, int j){
		Shape line = null;
		int gID, iNum, jNum, gIdNB, iNumNB, jNumNB; 
		float[] nbNode;
		ColorScale scale = new ColorScale();
		Color col = null;
		double xPos, yPos, xPosNB=0, yPosNB=0;
		double thetaRad, phiRad, thetaRadNB, phiRadNB;
//		double xPosINum = (this.voxelsize*translateXCoordRespective2Orig(Math.PI)) +this.border;
//		double yPosINum = (this.voxelsize*translateYCoordRespective2Orig(Math.PI/2)) +this.border;
//		double xPosINum = (translateXPixelCoordRespective2Orig(this.voxelsize*Math.PI)) +this.border;
//		double yPosINum = (translateYPixelCoordRespective2Orig(this.voxelsize*Math.PI/2)) +this.border;
		double xPosINum = getScreenPosFromRad(Math.PI, Math.PI/2).getFirst();
		double yPosINum = getScreenPosFromRad(Math.PI, Math.PI/2).getSecond();
		
		gID = (int) node[0];
		iNum = (int) node[1];
		jNum = (int) node[2];
		thetaRad = node[3];
		phiRad = node[4];				
		phiRad += Math.PI; // cause phi is expected to have a value of [0:pi]
		xPos = getScreenPosFromRad(phiRad, thetaRad).getFirst();
		yPos = getScreenPosFromRad(phiRad, thetaRad).getSecond();
		
		// --- gradient color edges between connected nodes
		nbNode = (float[]) this.nbhsNodes.get(j);
		thetaRadNB = nbNode[3];
		phiRadNB = nbNode[4];
		phiRadNB += Math.PI; // cause phi is expected to have a value of [0:pi]
		xPosNB = getScreenPosFromRad(phiRadNB, thetaRadNB).getFirst();
		yPosNB = getScreenPosFromRad(phiRadNB, thetaRadNB).getSecond();
		
		if (node[0]==nbNode[0] && node[1]==nbNode[1]){
//			System.out.print(nodeID +"\t" + lineID +"\t" + this.numNodesPerLine[lineID] +"\t");
			float ratio = (float)(nodeID+1)/(float)this.numNodesPerLine[lineID];
//			ratio = ((this.maxDistsJRes[lineID]-this.minDistsJRes[lineID])*node[2] + this.minDistsJRes[lineID])/this.maxDistsJRes[lineID];
//			System.out.print(ratio +"\t");
			
			// use ratio for color on scale
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
				int thres1 = 6; // 1-6:short range  6-25:middle range  25-n:long range
				int thres2 = 25;
				col = Color.black;
				if (Math.abs(jNum-iNum)<=thres1){
//					if (jNum-iNum < 0)
//						ratio = -1 * (float)Math.abs(jNum-iNum)/(float)(thres1);
//					else 
//						ratio = +1 * (float)Math.abs(jNum-iNum)/(float)(thres1);
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
					
					line = new Line2D.Double(xPos, yPos+this.yBorderThres, xPosINum, yPosINum+this.yBorderThres);		
					g2d.draw(line);
					line = new Line2D.Double(xPosINum, yPosINum+this.yBorderThres, xPosNB, yPosNB+this.yBorderThres);		
					g2d.draw(line);
				}
				else {
					// -- test for distance
//					if (Math.abs(xPosNB-xPos) > this.xDim/2){
//						if (xPos<xPosNB){
//							line = new Line2D.Double(xPos, yPos+this.yBorderThres, xPosNB-this.xDim, yPosNB+this.yBorderThres);		
//							g2d.draw(line);
//							line = new Line2D.Double(xPos+this.xDim, yPos+this.yBorderThres, xPosNB, yPosNB+this.yBorderThres);		
//							g2d.draw(line);
//						}
//						else{
//							line = new Line2D.Double(xPos-this.xDim, yPos+this.yBorderThres, xPosNB, yPosNB+this.yBorderThres);		
//							g2d.draw(line);
//							line = new Line2D.Double(xPos, yPos+this.yBorderThres, xPosNB+this.xDim, yPosNB+this.yBorderThres);		
//							g2d.draw(line);
//						}
//					}
//					else 
					{
						line = new Line2D.Double(xPos, yPos+this.yBorderThres, xPosNB, yPosNB+this.yBorderThres);		
						g2d.draw(line);
					}
				}	
			}									
			nodeID++;
		}		

//		// --- black edges between connected nodes
//		nbNode = (float[]) this.nbhsNodes.get(j);
//		xPosNB = (float) ((Math.PI+phiRadNB)*this.voxelsize);
//		yPosNB = thetaRadNB*this.voxelsize;
//		// if nodes have same gID and iNum --> connect through line
//		if (node[0]==nbNode[0] && node[1]==nbNode[1]){
//			g2d.setColor(Color.black);
//			line = new Line2D.Double(xPos, yPos, xPosNB, yPosNB);		
//			g2d.draw(line);
//		}
	}
	
	public void drawNBHSTraces(Graphics2D g2d){
		int gID, iNum, jResID;
		float[] node, nbNode;
		int lineID = 0;
		int nodeID = 0;

		int nbhsIndexC = 0;
		boolean specialRes = false;
				
		for(int i=0; i<this.nbhsNodes.size(); i++){
			
			node = (float[]) this.nbhsNodes.get(i);
//			System.out.println("0:graphi_id + '\t' + 1:i_num + '\t' + 2:j_num + '\t' + 3:theta + '\t' + 4:phi");
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
			if (this.nbhString.contains(String.valueOf(this.aas[jResID]))){
				specialRes = true;
			}
			// ordering is of importance
			if (nbhsIndexC > 0 && this.aas[jResID]==this.nbhsRes[nbhsIndexC-1])
				specialRes = true;
			if (this.aas[jResID] == this.nbhsRes[nbhsIndexC]){
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
				drawNBHSEdge(g2d, node, nodeID, lineID, j);
			}	
			else {
				lineID++;
//			    System.out.println();
				nodeID = 0;
			}		
		}	
		
	}
	
	protected void drawLongitudes(Graphics2D g2d){
		// Laengengrad
		g2d.setColor(this.longitudeColor);
		double xS, xE, yS, yE;
		Shape line;

		yS = ((0.0) *this.voxelsize) +this.border;
		yE = (Math.PI *this.voxelsize) +this.border;
//		yS = (translateYCoordRespective2Orig(0.0) *this.voxelsize) +this.border;
//		yE = (translateYCoordRespective2Orig(Math.PI) *this.voxelsize) +this.border;
		
		xS = (translateXCoordRespective2Orig(0.0) *this.voxelsize) +this.border;
		xE = xS;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
		xS = (translateXCoordRespective2Orig(Math.PI/2) *this.voxelsize) +this.border;
		xE = xS;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
		xS = (translateXCoordRespective2Orig(Math.PI) *this.voxelsize) +this.border;
		xE = xS;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
		xS = (translateXCoordRespective2Orig(3*Math.PI/2) *this.voxelsize) +this.border;
		xE = xS;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
	}
	
	protected void drawLatitudes(Graphics2D g2d){
		// Breitengrad
		g2d.setColor(this.latitudeColor);
		double xS, xE, yS, yE;
		Shape line;

		xS = ((0.0f) *this.voxelsize) +this.border;
		xE = (((2*Math.PI)) *this.voxelsize) +this.border;
//		xS = (translateYCoordRespective2Orig(0.0) *this.voxelsize) +this.border;
//		xE = (translateYCoordRespective2Orig(2*Math.PI) *this.voxelsize) +this.border;
		
		yS = (translateYCoordRespective2Orig(Math.PI/2) *this.voxelsize) +this.border;
		yE = yS;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
		yS = (translateYCoordRespective2Orig(Math.PI/4) *this.voxelsize) +this.border;
		yE = yS;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
		yS = (translateYCoordRespective2Orig(3*Math.PI/4) *this.voxelsize) +this.border;
		yE = yS;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
		
		yS = (translateYCoordRespective2Orig(Math.PI-0.01) *this.voxelsize) +this.border;
		yE = yS;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
		yS = (translateYCoordRespective2Orig(0+0.01) *this.voxelsize) +this.border;
		yE = yS;
		line = new Line2D.Double(xS,yS,xE,yE);
		g2d.draw(line);
	}
	
	protected void drawLongLatCentre(Graphics2D g2d){
		double bcenterx, bcentery;
		if (this.mapProjection){
			bcenterx = getScreenPosFromRad(this.origCoordinates.getFirst(), this.origCoordinates.getSecond()).getFirst();
			bcentery = getScreenPosFromRad(this.origCoordinates.getFirst(), this.origCoordinates.getSecond()).getSecond();
		}
		else {
			bcenterx = (this.origCoordinates.getFirst() * this.voxelsize) + this.border;
			bcentery = (this.origCoordinates.getSecond() * this.voxelsize) + this.border;
		}
		
		Shape circle = new Ellipse2D.Double(bcenterx-30, bcentery-30, 60, 60);
		g2d.draw(circle);
	}
	
	protected void drawCrosshair(Graphics2D g2d){
		// only in case of range selection we draw a diagonal cursor
		// drawing the cross-hair
		g2d.setColor(crosshairColor);
		g2d.drawLine(mousePos.x, 0, mousePos.x, g2dSize.height);
		g2d.drawLine(0, mousePos.y, g2dSize.width, mousePos.y);		
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
//		System.out.println(phiRanges.size()+"=?"+thetaRanges.size()+"=?"+selContacts.size());
//		System.out.println("Coordinates Occupied Rectangle:");
		Iterator<Pair<Float>> itrP = this.phiRanges.iterator();
		Iterator<Pair<Float>> itrT = this.thetaRanges.iterator();
//		Iterator<Pair<Integer>> itrC = selContacts.iterator();
		while (itrP.hasNext()){
			g2d.setColor(selAngleRangeColor);
			Pair<Float> phi = (Pair<Float>) itrP.next();
			Pair<Float> theta = (Pair<Float>) itrT.next();
//			Pair<Integer> con = (Pair<Integer>) itrC.next();
//			System.out.println("COR "+con.getFirst()+"_"+con.getSecond()+": "+phi.getFirst()+","+theta.getFirst()+" to "+phi.getSecond()+","+theta.getSecond());			
			float xS = phi.getFirst() * this.voxelsize;
			float xE = phi.getSecond() * this.voxelsize;
			float yS = theta.getFirst() * this.voxelsize;
			float yE = theta.getSecond() * this.voxelsize;
//			System.out.println(xS+","+yS+" to "+xE+","+yE);
			shape = new Rectangle2D.Float(xS, yS, xE-xS, yE-yS);
			g2d.draw(shape);
//			g2d.fill(shape);
		}
//		for(Pair<Float> phi:phiRanges) {}
		
	}
	
	private void drawSelectedAngleRange(){
		this.contStatBar.getAnglePanel().setIRes(String.valueOf(this.iRes));
		this.contStatBar.getAnglePanel().setJRes(String.valueOf(this.jRes));	
		
		float val = (float)(Math.round(this.tmpPhiRange.getFirst()*100))/100;
		this.contStatBar.getAnglePanel().setPhiMin(String.valueOf(val));	
		val = (float)(Math.round(this.tmpAngleComb.getFirst()*100))/100;
		this.contStatBar.getAnglePanel().setPhiMax(String.valueOf(val));	
		val = (float)(Math.round(this.tmpThetaRange.getFirst()*100))/100;
		this.contStatBar.getAnglePanel().setThetaMin(String.valueOf(val));	
		val = (float)(Math.round(this.tmpAngleComb.getSecond()*100))/100;
		this.contStatBar.getAnglePanel().setThetaMax(String.valueOf(val));
		
		this.contStatBar.repaint();		
	}
	
	// end drawing methods
	
	
	/**
	 * Returns the corresponding theta-phi values in the sphoxelView given screen
	 * coordinates
	 */
	private Pair<Float> screen2A(Point point){
		
		Pair<Float> floatP = new Pair<Float>((float)(point.x)/this.voxelsize,(float)(point.y)/this.voxelsize);
//		System.out.println("scree2A_point "+(point.x)+","+(point.y));
//		System.out.println("scree2A_pair "+(float)(point.x)/this.voxelsize+","+(float)(point.y)/this.voxelsize);
		
		return floatP;
	}
	
	public Pair<Double> translateCoordRespective2Orig(Pair<Double> pair){
		Pair<Double> transl;
		double xPos = translateXCoordRespective2Orig(pair.getFirst());
		double yPos = translateYCoordRespective2Orig(pair.getSecond());
		transl = new Pair<Double>(xPos, yPos);
		return transl;
	}
	
	public double translateXCoordRespective2Orig(double x){
		double dx = (this.origCoordinates.getFirst() - Math.PI);
		double xPos = x + dx;
//		xPos = (int)(xPos*1000)/1000;
		if (xPos > 2*Math.PI)
			xPos -= (2*Math.PI);
		else if (xPos < 0)
			xPos += (2*Math.PI);
		return xPos;
	}
	
	public double translateYCoordRespective2Orig(double y){
		double dy = (this.origCoordinates.getSecond() - (Math.PI/2));
		double yPos = y + dy;
//		yPos = (int)(yPos*1000)/1000;
//		if (yPos > Math.PI)
//			yPos -= Math.PI;
//		else if (yPos < 0)
//			yPos += Math.PI;
		return yPos;
	}
	
	public double translateXPixelCoordRespective2Orig(double x){
		double dx = (this.origCoordinates.getFirst() - Math.PI) * this.voxelsize;
		double xPos = x + dx;
		xPos = (int)(xPos*1000)/1000;
		if (xPos > 2*Math.PI* this.voxelsize)
			xPos -= (2*Math.PI* this.voxelsize);
		else if (xPos < 0)
			xPos += (2*Math.PI* this.voxelsize);
		return xPos;
	}
	
	public double translateYPixelCoordRespective2Orig(double y){
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
	 * Update tmpContact with the contacts contained in the rectangle given by
	 * the upperLeft and lowerRight.
	 */
	private void squareSelect(){
		Pair<Float> upperLeft = screen2A(mousePressedPos);
		Pair<Float> lowerRight = screen2A(mouseDraggingPos);
		// we reset the tmpContacts list so every new mouse selection starts
		// from a blank list
//		tmpContacts = new IntPairSet();

		float pmin = Math.min(upperLeft.getFirst(),lowerRight.getFirst());
		float tmin = Math.min(upperLeft.getSecond(),lowerRight.getSecond());
		float pmax = Math.max(upperLeft.getFirst(),lowerRight.getFirst());
		float tmax = Math.max(upperLeft.getSecond(),lowerRight.getSecond());

//		System.out.println("squareSelect "+pmin+"-"+pmax+" , "+tmin+"-"+tmax);
		
		this.tmpPhiRange = new Pair<Float>(pmin, pmax);
		this.tmpThetaRange = new Pair<Float>(tmin, tmax);
//		this.tmpSelContact = new Pair<Integer>(this.AAStr.indexOf(this.iRes),this.AAStr.indexOf(this.jRes));
		this.tmpSelContact = new Pair<Integer>(this.iNum,this.jNum);
		
		this.tmpAngleComb = new Pair<Float>(pmax, tmax);

	}
	
//	/** Called by ResidueRuler to enable display of ruler "crosshair" */	
//	public void showRulerCrosshair() {
//		showRulerCrosshair = true;
//	}
//	/** Called by ResidueRuler to switch off showing ruler coordinates */
//	public void hideRulerCoordinate() {
//		showRulerCoord = false;
//	}

	/**
	 * Main method to draw the component on screen. This method is called each
	 * time the component has to be (re) drawn on screen. It is called
	 * automatically by Swing or by explicitly calling cmpane.repaint().
	 */
	@Override
	protected synchronized void paintComponent(Graphics g) {
		Graphics2D g2d = (Graphics2D) g.create();
		
		this.screenSize = this.contactView.getScreenSize();
//		System.out.println("ContactPane paintComponent HxW: "+this.screenSize.height+"x"+this.screenSize.width);
//		System.out.println("ContactPane pC HxW: "+this.getSize().height+"x"+this.getSize().width+"  --> ratio="+(float)(this.getSize().height)/(float)(this.getSize().width));
		
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
		
//		System.out.println("CPane g2dSize HxB: "+this.g2dSize.height+"x"+this.g2dSize.width);
		
//		Dimension dim = this.contactView.getSize();		
//		this.xDim = (int) dim.getWidth();
//		this.yDim = (int) dim.getHeight();
		
		// update pixel dimensions
		this.pixelWidth = (float)(this.xDim-2*this.border)/(float)(2*this.numSteps) ;
		this.pixelHeight =  this.pixelWidth; //(this.yDim-2*this.border-this.yBorderThres)/(2*this.numSteps);
		this.voxelsize = (float) ((float)(this.numSteps*this.pixelWidth)/Math.PI);
//		System.out.println("NumSteps= "+this.numSteps+"  PixelSize: "+this.pixelHeight+"x"+this.pixelWidth+" VoxelS:"+this.voxelsize);
		
		g2d.setBackground(this.backgroundColor);
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
//		drawSphoxels(g2d);
		drawSphoxelMap(g2d);
		drawNBHSTraces(g2d);
		// drawing selection rectangle if dragging mouse and showing temp
		// selection in red (tmpContacts)
		if (dragging && contactView.getGUIState().getSelectionMode()==ContactGUIState.SelMode.RECT) {
			g2d.setColor(squareSelColor);
			int xmin = Math.min(mousePressedPos.x,mouseDraggingPos.x);
			int ymin = Math.min(mousePressedPos.y,mouseDraggingPos.y);
			int xmax = Math.max(mousePressedPos.x,mouseDraggingPos.x);
			int ymax = Math.max(mousePressedPos.y,mouseDraggingPos.y);
			g2d.drawRect(xmin,ymin,xmax-xmin,ymax-ymin);
		} 
		drawCrosshair(g2d);
		drawLongitudes(g2d);
		drawLatitudes(g2d);
		drawLongLatCentre(g2d);
		drawOccupiedAngleRanges(g2d);
		drawSelectedAngleRange();
	}
	
	/*---------------------------- setters and getters -----------------------------*/
	
	public void setStatusBar(ContactStatusBar statusBar) {
		this.contStatBar = statusBar;		
	}
	public ContactStatusBar getSatusBar(){
		return this.contStatBar;
	}
	
	public Dimension getScreenSize(){
		return this.screenSize;
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
	}
	public void setResol(float res){
		this.numSteps = (int) (180/res);
		this.resol = 180/(float)this.numSteps;
		System.out.println("numSteps="+this.numSteps+"  --> resol="+this.resol);
	}	
	public Pair<Float> getOrigCoord(){
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
	}

	public double getMaxAllowedRat() {
		return maxAllowedRat;
	}

	public void setMaxAllowedRat(double maxAllowedRat) {
		this.maxAllowedRat = maxAllowedRat;
	}


	public int getChosenColourScale() {
		return chosenColourScale;
	}

	public void setChosenColourScale(int chosenColourScale) {
		this.chosenColourScale = chosenColourScale;
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
	
	/*---------------------------- mouse events -----------------------------*/

	public void mousePressed(MouseEvent evt) {
		// This is called when the user presses the mouse anywhere
		// in the frame
		
		System.out.println("mousePressed");

		lastMouseButtonPressed = evt.getButton();
		mousePressedPos = evt.getPoint();
		if(lastMouseButtonPressed == MouseEvent.BUTTON2) dragging = false;
	}

	public void mouseReleased(MouseEvent evt) {

		// Called whenever the user releases the mouse button.
		// TODO: Move much of this to MouseClicked and pull up Contact cont = screen2A...
		if (evt.isPopupTrigger()) {
			dragging = false;
			return;
		}
		// only if release after left click (BUTTON1)
		if (evt.getButton()==MouseEvent.BUTTON1) {

			switch (contactView.getGUIState().getSelectionMode()) {
			case RECT:
				if (dragging){
//					squareSelect();
					if (this.selContacts.size()>0 && 
							this.selContacts.lastElement().getFirst()==this.tmpSelContact.getFirst() && 
							this.selContacts.lastElement().getSecond()==this.tmpSelContact.getSecond()){
						this.phiRanges.removeElementAt(this.phiRanges.size()-1);
						this.thetaRanges.removeElementAt(this.thetaRanges.size()-1);
						this.selContacts.removeElementAt(this.selContacts.size()-1);						
					}
					this.phiRanges.add(this.tmpPhiRange);
					this.thetaRanges.add(this.tmpThetaRange);
					this.selContacts.add(this.tmpSelContact);
					updateAngleRange();
					System.out.println("mouseReleased "+this.tmpPhiRange.getFirst()+"-"+this.tmpPhiRange.getSecond()
							+" , "+this.tmpThetaRange.getFirst()+"-"+this.tmpThetaRange.getSecond());
				}				
				dragging = false;
				this.repaint();
				return;
				
			case CLUSTER:
//				// resets selContacts when clicking mouse

//				this.repaint();
				return;
				
			case PAN:
				float resol = (float) (Math.PI/this.numSteps);
				int fac = (int) Math.ceil(this.tmpAngleComb.getFirst()/resol);
				float first = fac*resol; //this.tmpAngleComb.getFirst();
				fac = (int) Math.ceil(this.tmpAngleComb.getSecond()/resol);
				float second = fac*resol;
				
				second = (float) (Math.PI/2);   // --> panning just horizontally
				this.origCoordinates = new Pair<Float>(first,second); // this.tmpAngleComb;
				this.contactView.phiRuler.repaint();
				this.contactView.thetaRuler.repaint();
				System.out.println("tmpAngleComb: "+this.tmpAngleComb.getFirst()+","+this.tmpAngleComb.getSecond());
				System.out.println("OrigCoord: "+this.origCoordinates.getFirst()+","+this.origCoordinates.getSecond());
				this.repaint();
				return;
			}
			
			this.tmpPhiRange = new Pair<Float>(0.0f, 0.0f);
			this.tmpThetaRange = new Pair<Float>(0.0f, 0.0f);
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
				squareSelect();
				break;
			case PAN:
				float resol = (float) (Math.PI/this.numSteps);
				int fac = (int) Math.ceil(this.tmpAngleComb.getFirst()/resol);
				float first = fac*resol; //this.tmpAngleComb.getFirst();
				fac = (int) Math.ceil(this.tmpAngleComb.getSecond()/resol);
				float second = fac*resol;
				
				second = (float) (Math.PI/2);   // --> panning just horizontally
				this.origCoordinates = new Pair<Float>(first,second); // this.tmpAngleComb;
				this.contactView.phiRuler.repaint();
				this.contactView.thetaRuler.repaint();
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
		this.repaint();
	}

	public void mouseClicked(MouseEvent evt) {
	}

	public void mouseMoved(MouseEvent evt) {
//		System.out.println("mouseMoved");

		mousePos = evt.getPoint();
		
		Pair<Float> pos = screen2A(mousePos);
		this.tmpAngleComb = new Pair<Float>(pos.getFirst(), pos.getSecond());
		this.repaint();
	}

	public void keyPressed(KeyEvent e) {
		// TODO Auto-generated method stub
		System.out.println("keyPressed");
		//-- Process arrow "virtual" keys
    	float first = this.origCoordinates.getFirst();
    	float second = this.origCoordinates.getSecond();
    	float fac = (float) (Math.PI / 20);
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
        	first += Math.PI;
        if (second>Math.PI)
        	first -= Math.PI;
        
        second = (float) (Math.PI/2);   // --> panning just horizontally     	
		this.origCoordinates = new Pair<Float>(first,second); // this.tmpAngleComb;
		this.contactView.phiRuler.repaint();
		this.contactView.thetaRuler.repaint();		
	}

	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub
		System.out.println("keyReleased");
		int id = e.getID();
        String keyString;
        if (id == KeyEvent.KEY_TYPED) {
            char c = e.getKeyChar();
            keyString = "key character = '" + c + "'";
        } else {
            int keyCode = e.getKeyCode();
            keyString = "key code = " + keyCode
                    + " ("
                    + KeyEvent.getKeyText(keyCode)
                    + ")";
        }
        System.out.println("keyString= "+keyString);
	}

	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub	
		System.out.println("keyReleased");	
	}

}
