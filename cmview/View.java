package cmview;
import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.awt.image.BufferedImage;
import javax.imageio.*;
import java.util.HashMap;

import cmview.datasources.ContactMapFileModel;
import cmview.datasources.GraphDbModel;
import cmview.datasources.Model;
import cmview.datasources.MsdsdModel;
import cmview.datasources.PdbFileModel;
import cmview.datasources.PdbaseModel;
import proteinstructure.Contact;

/**
 * The main GUI window.
 * 
 * Initialized with mod=null, an empty window with the menu bars is shown.
 * Initialized with a valid model, the contact map is displayed in a ContactMapPane.
 * Multiple instances of this class are created to show multiple contact maps.
 * 
 * @author	Juliane Dinse
 * @author	Henning Stehr
 * @author	Jose Duarte
 * 
 * Class: 	View
 * Package:	cmview
 * Date:	20/02/2007 last update 12/06/2007
 *
 */
public class View extends JFrame implements ActionListener {

	static final long serialVersionUID = 1l;
	protected static final int SQUARE_SEL = 1;
	protected static final int FILL_SEL = 2;
	protected static final int NODE_NBH_SEL = 3;
	protected static final int SHOW_COMMON_NBH = 4;
	
	
	// GUI components
	JLabel statusPane; 			// status bar
	JPanel cmp; 				// contact Map Panel
	JPanel topRul;
	JPanel leftRul;
	JPopupMenu popup; 		// context menu
	JFileChooser fileChooser;
	JColorChooser colorChooser;

	// Menu items
	JMenuItem sendM, squareM, fillM, loadPDBM, comNeiM, triangleM, nodeNbhSelM;
	JMenuItem sendP, squareP, fillP, loadPDBP, comNeiP, triangleP, nodeNbhSelP;
	JMenuItem mmLoadGraph, mmLoadPdbase, mmLoadMsd, mmLoadCm, mmLoadPdb;
	JMenuItem mmSaveGraph, mmSaveCm, mmSavePng;
	JMenuItem mmViewShowPdbResSers, mmViewHighlightComNbh, mmViewRulers;
	JMenuItem mmSelectAll;
	JMenuItem mmColorReset, mmColorPaint, mmColorChoose;
	JMenuItem mmInfo, mmPrint, mmQuit, mmHelpAbout, mmHelpHelp;


	// Data and status variables

	private Model mod;
	public ContactMapPane cmPane;
	public ResidueRuler topRuler;
	public ResidueRuler leftRuler;
	private PyMolAdaptor pymolAdaptor;
	private String pyMolServerUrl;
	private int currentAction;
	private int pymolSelSerial;		 	// for incrementation numbering TODO: move to Model
	private int pymolNbhSerial;

	private boolean doShowPdbSers;
	private boolean highlightComNbh;
	private boolean showRulers;
	private Color currentPaintingColor;

	private HashMap<Contact,Integer> comNbhSizes;

	/** Create a new View object */
	public View(Model mod, String title, String pyMolServerUrl) {
		super(title);
		this.mod = mod;
		this.pyMolServerUrl=pyMolServerUrl;
		if(mod == null) {
			this.setPreferredSize(new Dimension(800,800));
		}
		this.initGUI();
		this.currentAction = SQUARE_SEL;
		this.pymolSelSerial = 1;
		this.pymolNbhSerial = 1;
		this.doShowPdbSers = false;
		this.highlightComNbh = false;
		this.showRulers=false;
		this.currentPaintingColor = Color.blue;
		
		fileChooser = new JFileChooser();
		colorChooser = new JColorChooser();
		colorChooser.setPreviewPanel(new JPanel()); // removing the preview panel
	}

	/** Initialize and show the main GUI window */
	private void initGUI(){

		// Setting the main layout 
		setLayout(new BorderLayout());
		//setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setLocation(20,20);

		// Creating the Panels
		statusPane = new JLabel("Click right mouse button for context menu");
		cmp = new JPanel(new BorderLayout()); // Contact Map Panel
		topRul = new JPanel(new BorderLayout());
		leftRul = new JPanel(new BorderLayout());
		
		// Adding the context menu
		JPopupMenu.setDefaultLightWeightPopupEnabled(false);
		popup = new JPopupMenu();

		ImageIcon icon1 = new ImageIcon("icons/shape_square.png");
		ImageIcon icon2 = new ImageIcon("icons/paintcan.png");
		ImageIcon icon3 = new ImageIcon("icons/group.png");
		ImageIcon icon4 = new ImageIcon("icons/shape_square_go.png");
		ImageIcon icon5 = new ImageIcon("icons/shape_flip_horizontal.png");
		ImageIcon icon6 = new ImageIcon("icons/shape_rotate_clockwise.png");
		//ImageIcon icon7 = new ImageIcon("icons/picture_go.png");

//		ImageIcon icon8 = new ImageIcon("icons/information.png");
//		ImageIcon icon9 = new ImageIcon("icons/printer.png");
//		ImageIcon icon10 = new ImageIcon("icons/door_open.png");

		squareP = new JMenuItem("Square Selection Mode", icon1);
		fillP = new JMenuItem("Fill Selection Mode", icon2);
		nodeNbhSelP = new JMenuItem("Node Neighbourhood Selection Mode", icon3);
		sendP = new JMenuItem("Show Selected Edges in PyMol", icon4);
		comNeiP = new JMenuItem("Show Common Neighbours Mode", icon5);
		triangleP = new JMenuItem("Show Common Neighbour Triangles in PyMol", icon6);

		squareP.addActionListener(this);
		fillP.addActionListener(this);
		nodeNbhSelP.addActionListener(this);
		comNeiP.addActionListener(this);
		sendP.addActionListener(this);
		triangleP.addActionListener(this);		

		popup.add(squareP);
		popup.add(fillP);
		popup.add(nodeNbhSelP);
		popup.addSeparator();	
		popup.add(sendP);
		popup.addSeparator();	
		popup.add(comNeiP);
		popup.add(triangleP);

		if(mod != null) {
			cmPane = new ContactMapPane(mod, this);			
			cmp.add(cmPane);
			topRuler = new ResidueRuler(cmPane,mod,this,ResidueRuler.TOP);
			leftRuler = new ResidueRuler(cmPane,mod,this,ResidueRuler.LEFT);
			topRul.add(topRuler);
			leftRul.add(leftRuler);
		}

		// Creating the menu bar
		JMenuBar menuBar;
		JMenu menu, submenu;

		menuBar = new JMenuBar();

		// File menu
		menu = new JMenu("File");
		menu.setMnemonic(KeyEvent.VK_F);
		mmInfo = new JMenuItem("Info");
		mmInfo.addActionListener(this);
		menu.add(mmInfo);
		submenu = new JMenu("Load from");
		submenu.setMnemonic(KeyEvent.VK_L);

		mmLoadGraph = new JMenuItem("Graph database");
		mmLoadPdbase = new JMenuItem("Pdbase");
		mmLoadMsd = new JMenuItem("MSD");
		mmLoadPdb = new JMenuItem("PDB file");
		mmLoadCm = new JMenuItem("Contact map file");

		submenu.add(mmLoadGraph);
		submenu.add(mmLoadPdbase);
		submenu.add(mmLoadMsd);
		submenu.add(mmLoadPdb);
		submenu.add(mmLoadCm);

		mmLoadGraph.addActionListener(this);
		mmLoadPdbase.addActionListener(this);			
		mmLoadMsd.addActionListener(this);			
		mmLoadPdb.addActionListener(this);			
		mmLoadCm.addActionListener(this);			

		menu.add(submenu);
		submenu = new JMenu("Save to");
		submenu.setMnemonic(KeyEvent.VK_S);

		mmSaveGraph = new JMenuItem("Graph database");			
		mmSaveCm = new JMenuItem("Contact map file");
		mmSavePng = new JMenuItem("PNG file");

		//submenu.add(mmSaveGraph);
		submenu.add(mmSaveCm);			
		submenu.add(mmSavePng);

		mmSaveGraph.addActionListener(this);			
		mmSaveCm.addActionListener(this);			
		mmSavePng.addActionListener(this);	

		menu.add(submenu);
		mmPrint = new JMenuItem("Print");
		mmPrint.addActionListener(this);
		menu.add(mmPrint);
		mmQuit = new JMenuItem("Quit");
		mmQuit.addActionListener(this);
		menu.add(mmQuit);
		menuBar.add(menu);

		// View menu
		menu = new JMenu("View");
		menu.setMnemonic(KeyEvent.VK_V);
		mmViewShowPdbResSers = new JMenuItem("Toggle show PDB residue numbers");
		mmViewHighlightComNbh = new JMenuItem("Toggle highlight of cells by common neighbourhood size");
		mmViewRulers = new JMenuItem("Toggle rulers");
		menu.add(mmViewShowPdbResSers);
		menu.add(mmViewRulers);
		menu.addSeparator();
		menu.add(mmViewHighlightComNbh);
		mmViewShowPdbResSers.addActionListener(this);
		mmViewHighlightComNbh.addActionListener(this);
		mmViewRulers.addActionListener(this);
		menuBar.add(menu);

		// Select menu
		menu = new JMenu("Select");
		menu.setMnemonic(KeyEvent.VK_S);
		mmSelectAll = new JMenuItem("All contacts");
		menu.add(mmSelectAll);
		mmSelectAll.addActionListener(this);
		menuBar.add(menu);
		
		// Color menu
		menu = new JMenu("Color");
		menu.setMnemonic(KeyEvent.VK_C);
		mmColorReset= new JMenuItem("Reset contact colors to black");
		mmColorChoose = new JMenuItem("Choose painting color");
		mmColorPaint = new JMenuItem("Paint current selection");
		menu.add(mmColorChoose);
		menu.add(mmColorPaint);
		menu.add(mmColorReset);
		mmColorReset.addActionListener(this);
		mmColorPaint.addActionListener(this);
		mmColorChoose.addActionListener(this);
		menuBar.add(menu);
		
		// Action menu
		menu = new JMenu("Action");
		menu.setMnemonic(KeyEvent.VK_A);

		squareM = new JMenuItem("Square Selection Mode", icon1);
		fillM = new JMenuItem("Fill Selection Mode", icon2);
		nodeNbhSelM = new JMenuItem("Node Neighbourhood Selection Mode", icon3);
		sendM = new JMenuItem("Show Selected Edges in PyMol", icon4);
		comNeiM = new JMenuItem("Show Common Neighbours Mode", icon5);
		triangleM = new JMenuItem("Show Common Neighbour Triangles in PyMol", icon6);

		squareM.addActionListener(this);
		fillM.addActionListener(this);
		nodeNbhSelM.addActionListener(this);
		comNeiM.addActionListener(this);
		sendM.addActionListener(this);
		triangleM.addActionListener(this);			

		menu.add(squareM);
		menu.add(fillM);
		menu.add(nodeNbhSelM);
		menu.addSeparator();
		menu.add(sendM);
		menu.addSeparator();		
		menu.add(comNeiM);
		menu.add(triangleM);

		menuBar.add(menu);

		// Help menu
		menu = new JMenu("Help");
		menu.setMnemonic(KeyEvent.VK_H);	
		mmHelpAbout = new JMenuItem("About");
		mmHelpHelp = new JMenuItem("Help");
		mmHelpAbout.addActionListener(this);
		mmHelpHelp.addActionListener(this);
		menu.add(mmHelpHelp);
		menu.add(mmHelpAbout);
		menuBar.add(menu);

		this.setJMenuBar(menuBar);
		this.getContentPane().add(cmp,BorderLayout.CENTER);
		if(showRulers) {
			this.getContentPane().add(topRul, BorderLayout.NORTH);
			this.getContentPane().add(leftRul, BorderLayout.WEST);
		}
		//this.getContentPane().add(statusPane,BorderLayout.SOUTH);

		// Show GUI
		pack();
		setVisible(true);

	}

	/*------------------------------ event handling -------------------------*/

	public void actionPerformed (ActionEvent e) {

		// Action menu

		// square button clicked
		if (e.getSource() == squareM || e.getSource() == squareP) {

			currentAction = SQUARE_SEL;

		}
		// fill button clicked
		if (e.getSource() == fillM || e.getSource() == fillP) {

			currentAction = FILL_SEL;
		}
		// node neihbourhood selection button clicked 
		if (e.getSource() == nodeNbhSelM || e.getSource() == nodeNbhSelP ) {

			currentAction = NODE_NBH_SEL;
		}		
		// showing com. Nei. button clicked
		if (e.getSource() == comNeiM || e.getSource() == comNeiP) {

			currentAction = SHOW_COMMON_NBH;
		}
		// send selection button clicked
		if (e.getSource() == sendM || e.getSource() == sendP) {	
			if(mod==null) {
				showNoContactMapWarning();
			} else {
				pymolAdaptor.edgeSelection(pymolSelSerial, cmPane.getSelContacts());
				this.pymolSelSerial++;
			}
		}
		// send com.Nei. button clicked
		if(e.getSource()== triangleM || e.getSource()== triangleP) {
			if(mod==null) {
				showNoContactMapWarning();
			} else {
				pymolAdaptor.showTriangles(cmPane.getCommonNbh(),pymolNbhSerial);
				this.pymolNbhSerial++;
			}
		}

		// File Menu
		// Load
		if(e.getSource() == mmLoadGraph) {
			handleLoadFromGraphDb();
		}
		if(e.getSource() == mmLoadPdbase) {
			handleLoadFromPdbase();
		}		  
		if(e.getSource() == mmLoadMsd) {
			handleLoadFromMsd();
		}			  
		if(e.getSource() == mmLoadPdb) {
			handleLoadFromPdbFile();
		}
		if(e.getSource() == mmLoadCm) {
			handleLoadFromCmFile();
		}
		// Save
		if(e.getSource() == mmSaveGraph) {
			handleSaveToGraphDb();
		}		  
		if(e.getSource() == mmSaveCm) {
			handleSaveToCmFile();
		}
		if(e.getSource() == mmSavePng) {
			handleSaveToPng();
		}	
		// Info, Print, Quit
		if(e.getSource() == mmInfo) {
			handleInfo();
		}		  		  
		if(e.getSource() == mmPrint) {
			handlePrint();
		}		  
		if(e.getSource() == mmQuit) {
			handleQuit();
		}

		// View Menu
		if(e.getSource() == mmColorReset) {
			handleViewReset();
		}
		if(e.getSource() == mmColorPaint) {
			handleViewColor();
		}
		if(e.getSource() == mmColorChoose) {
			handleViewSelectColor();
		}
		if(e.getSource() == mmViewShowPdbResSers) {
			doShowPdbSers = !doShowPdbSers;
		}		  
		if(e.getSource() == mmViewHighlightComNbh) {
			if(mod==null) {
				showNoContactMapWarning();
			} else {
				highlightComNbh = !highlightComNbh;
				if (highlightComNbh) comNbhSizes = mod.getAllEdgeNbhSizes();
				cmPane.repaint();
			}
		}		  
		if(e.getSource() == mmViewRulers) {
			toggleRulers();
		}		  
		
		// Select menu
		if(e.getSource() == mmSelectAll) {
			handleSelectAll();
		}		
		
		// Help Menu
		if(e.getSource() == mmHelpAbout) {
			handleHelpAbout();
		}
		if(e.getSource() == mmHelpHelp) {
			handleHelpHelp();
		}		

	}

	private void handleLoadFromGraphDb() {
		LoadDialog dialog = new LoadDialog(this, "Load from graph database", new LoadAction() {
			public void doit(Object o, String f, String ac, String cc, String ct, double dist, int minss, int maxss, String db, int gid) {
				View view = (View) o;
				view.doLoadFromGraphDb(db, gid);
			}
		}, null, null, null, null, null, null, null, "pdb_reps_graph", "");
		dialog.showIt();
	}

	public void doLoadFromGraphDb(String db, int gid) {
		System.out.println("Loading from graph database");
		System.out.println("Database:\t" + db);
		System.out.println("Graph Id:\t" + gid);
		Model mod = new GraphDbModel(gid, db);
		this.spawnNewViewWindow(mod);
	}

	private void handleLoadFromPdbase() {
//		String pdbCode = "1tdr";
//		String chainCode = "B";
//		String edgeType = "Ca";
//		double distCutoff = 8.0;
//		int seqSep = 0;
		LoadDialog dialog = new LoadDialog(this, "Load from Pdbase", new LoadAction() {
			public void doit(Object o, String f, String ac, String cc, String ct, double dist, int minss, int maxss, String db, int gid) {
				View view = (View) o;
				view.doLoadFromPdbase(ac, cc, ct, dist, minss, maxss, db);
			}
		}, null, "", "", "", "", "", "", null, null);
		dialog.showIt();		  

	}

	public void doLoadFromPdbase(String ac, String cc, String ct, double dist, int minss, int maxss, String db) {
		System.out.println("Loading from Pdbase");
		System.out.println("PDB code:\t" + ac);
		System.out.println("Chain code:\t" + cc);
		System.out.println("Contact type:\t" + ct);
		System.out.println("Dist. cutoff:\t" + dist);	
		System.out.println("Min. Seq. Sep.:\t" + minss);
		System.out.println("Max. Seq. Sep.:\t" + maxss);
		db = "pdbase";
		System.out.println("Database:\t" + db);	
		Model mod = new PdbaseModel(ac, cc, ct, dist, minss, maxss, db);
		this.spawnNewViewWindow(mod);			  
	}

	private void handleLoadFromMsd() {
//		String pdbCode = "1tdr";
//		String chainCode = "B";
//		String edgeType = "Ca";
//		double distCutoff = 8.0;
//		int seqSep = 0;
		LoadDialog dialog = new LoadDialog(this, "Load from MSD", new LoadAction() {
			public void doit(Object o, String f, String ac, String cc, String ct, double dist, int minss, int maxss, String db, int gid) {
				View view = (View) o;
				view.doLoadFromMsd(ac, cc, ct, dist, minss, maxss, db);
			}
		}, null, "", "", "", "", "", "", null, null);
		dialog.showIt();
	}

	public void doLoadFromMsd(String ac, String cc, String ct, double dist, int minss, int maxss, String db) {
		System.out.println("Loading from MSD");
		System.out.println("PDB code:\t" + ac);
		System.out.println("Chain code:\t" + cc);
		System.out.println("Contact type:\t" + ct);
		System.out.println("Dist. cutoff:\t" + dist);	
		System.out.println("Min. Seq. Sep.:\t" + minss);
		System.out.println("Max. Seq. Sep.:\t" + maxss);
		db = "msdsd_00_07_a";
		System.out.println("Database:\t" + db);	
		Model mod = new MsdsdModel(ac, cc, ct, dist, minss, maxss, db);
		this.spawnNewViewWindow(mod);			  
	}	  

	private void handleLoadFromPdbFile() {
//		String chainCode = "A";
//		String edgeType = "Ca";
//		double distCutoff = 8.0;
//		int seqSep = 0;
		LoadDialog dialog = new LoadDialog(this, "Load from Pdb file", new LoadAction() {
			public void doit(Object o, String f, String ac, String cc, String ct, double dist, int minss, int maxss, String db, int gid) {
				View view = (View) o;
				view.doLoadFromPdbFile(f, cc, ct, dist, minss, maxss);
			}
		}, "", null, "", "", "", "", "", null, null);
		dialog.showIt();
	}

	public void doLoadFromPdbFile(String f, String cc, String ct, double dist, int minss, int maxss) {
		System.out.println("Loading from Pdb file");
		System.out.println("Filename:\t" + f);
		System.out.println("Chain code:\t" + cc);
		System.out.println("Contact type:\t" + ct);
		System.out.println("Dist. cutoff:\t" + dist);	
		System.out.println("Min. Seq. Sep.:\t" + minss);
		System.out.println("Max. Seq. Sep.:\t" + maxss);
		Model mod = new PdbFileModel(f, cc, ct, dist, minss, maxss);
		this.spawnNewViewWindow(mod);
	}	

	private void handleLoadFromCmFile() {
		LoadDialog dialog = new LoadDialog(this, "Load from Contact map file", new LoadAction() {
			public void doit(Object o, String f, String ac, String cc, String ct, double dist, int minss, int maxss, String db, int gid) {
				View view = (View) o;
				view.doLoadFromCmFile(f);
			}
		}, "", null, null, null, null, null, null, null, null);
		dialog.showIt();		  
	}

	public void doLoadFromCmFile(String f) {
		System.out.println("Loading from Pdb file");
		System.out.println("Filename:\t" + f);
		Model mod = new ContactMapFileModel(f);
		this.spawnNewViewWindow(mod);
	}		  

	private void handleSaveToGraphDb() {
		System.out.println("Saving to graph db not implemented yet");
	}

	public void doSaveToGraphDb() {
		// TODO
	}

	private void handleSaveToCmFile() {
		if(this.mod == null) {
			//System.out.println("No contact map loaded yet.");
			showNoContactMapWarning();
		} else {
			int ret = fileChooser.showSaveDialog(this);
			if(ret == JFileChooser.APPROVE_OPTION) {
				File chosenFile = fileChooser.getSelectedFile();
				String path = chosenFile.getPath();
				try {
					this.mod.writeToContactMapFile(path);
				} catch(IOException e) {
					System.err.println("Error writing to file " + path);
				}
			}
		}
	}

	private void handleSaveToPng() {
		if(this.mod == null) {
			//System.out.println("No contact map loaded yet.");
			showNoContactMapWarning();
		} else {
			int ret = fileChooser.showSaveDialog(this);
			if(ret == JFileChooser.APPROVE_OPTION) {
				File chosenFile = fileChooser.getSelectedFile();

				// Create a buffered image in which to draw
				BufferedImage bufferedImage = new BufferedImage(cmPane.getWidth(), cmPane.getHeight(), BufferedImage.TYPE_INT_RGB);

				// Create a graphics contents on the buffered image
				Graphics2D g2d = bufferedImage.createGraphics();

				// Draw the current contact map window to Image
				cmPane.paintComponent(g2d);

				try {
					ImageIO.write(bufferedImage, "png", chosenFile);
					System.out.println("File " + chosenFile.getPath() + " saved.");
				} catch (IOException e) {
					System.err.println("Error while trying to write to PNG file " + chosenFile.getPath());
				}
			}
		}
	}

	private void handleInfo() {
		if(this.mod == null) {
			//System.out.println("No contact map loaded yet.");
			showNoContactMapWarning();
		} else {
			String seq = mod.getSequence();
			String s = seq.length() <= 10?(seq.length()==0?"Unknown":seq):seq.substring(0,10) + "...";
			String message = "Pdb code: " + mod.getPDBCode() + "\n"
			+ "Chain code: " + mod.getChainCode() + "\n"
			+ "Contact type: " + mod.getContactType() + "\n"
			+ "Distance cutoff: " + mod.getDistanceCutoff() + "\n"
			+ "Min Seq Sep: " + mod.getMinSequenceSeparation() + "\n"
			+ "Max Seq Sep: " + mod.getMaxSequenceSeparation() + "\n"
			+ "\n"
			+ "Contact map size: " + mod.getMatrixSize() + "\n"
			+ "Unobserved Residues: " + mod.getNumberOfUnobservedResidues() + "\n"
			+ "Number of contacts: " + mod.getNumberOfContacts() + "\n"
			+ "Directed: " + (mod.isDirected()?"Yes":"No")
			+ "\n"
			+ "Sequence: " + s;
			JOptionPane.showMessageDialog(this,
					message,
					"Contact map info",
					JOptionPane.PLAIN_MESSAGE);
		}
	}

	private void handlePrint() {
		if(this.mod == null) {
			//System.out.println("No contact map loaded yet.");
			showNoContactMapWarning();
		} else {
			PrintUtil.printComponent(this.cmPane);
		}
	}

	private void handleQuit() {
		System.exit(0);
	}

	private void handleViewReset() {
		//System.out.println("View:Reset not implemented yet");
		if(mod==null) {
			showNoContactMapWarning();
		} else {
			cmPane.resetColorMap();
		}
	}	  

	private void handleViewColor() {
		//System.out.println("View:Color by type not implemented yet");
		if(mod==null) {
			showNoContactMapWarning();
		} else {
			cmPane.paintCurrentSelection(currentPaintingColor);
			cmPane.resetSelContacts();
			cmPane.repaint();
		}
	}
	
	private void handleViewSelectColor() {	
		JDialog chooseDialog = JColorChooser.createDialog(this, "Choose color", true, colorChooser, 
				new ActionListener() {
					public void actionPerformed(ActionEvent e){
							Color c = colorChooser.getColor();
							currentPaintingColor = c;
							System.out.println("Color changed to " + c);
					}
				} , null); // no handler for cancel button
		chooseDialog.setVisible(true);
	}	

	private void handleSelectAll() {
		cmPane.selectAllContacts();
	}
	
	private void handleHelpAbout() {
		JOptionPane.showMessageDialog(this,
				"           Contact Map Viewer v"+Start.VERSION+"\n" + 
				"               (C) AG Lappe 2007",
				"About",
				JOptionPane.PLAIN_MESSAGE);
	}
	
	private void handleHelpHelp() {
		JOptionPane.showMessageDialog(this,
				"General\n" +
				"- Click right mouse button in contact map for a context menu of available actions\n" +
				"\n" +
				"Square selection mode\n" +
				"- Click on a contact to select it\n" +
				"- Drag the mouse to select a rectangular area of contacts\n" +
				"- Hold 'Ctrl' while selecting to add to the current selection\n" +
				"- Click on a non-contact to reset the current selection\n" +
				"\n" +
				"Fill selection mode\n" +
				"- Click on a contact to start a fill selection from that contact\n" +
				"- Hold 'Ctrl' while selecting to add to the current selection\n" +
				"\n" +
				"Node neighbourhood selection mode\n" +
				"- Click on a residue in the ruler or in the diagonal to select its contacts\n" +
				"- Click on a cell in the upper half to select all contacts of that pair of residues\n" +
				"- Hold 'Ctrl' while selecting to add to the current selection\n" +				
				"\n" +
				"Show selected contacts in pymol\n" +
				"- Shows the currently selected contacts as edges in Pymol\n" +
				"\n" +
				"Show common neigbours\n" +
				"- Click on a contact or non-contact to see the common neighbours for that pair of residues\n" +
				"\n" +
				"Show common neighbours in pymol\n" +
				"- Shows the last shown common neighbours as triangles in pymol\n",
				"Help",
				JOptionPane.PLAIN_MESSAGE);
	}	
	
	/** Shows a window with a warning message that no contact map is loaded yet */
	private void showNoContactMapWarning() {
		JOptionPane.showMessageDialog(this, "No contact map loaded yet", "Error", JOptionPane.INFORMATION_MESSAGE);
		
	}

	/*---------------------------- public methods ---------------------------*/

	/** show/hide rulers and toggle the value of showRulers */
	private void toggleRulers() {
		showRulers = !showRulers;
		if(showRulers) {
			this.getContentPane().add(topRul, BorderLayout.NORTH);
			this.getContentPane().add(leftRul, BorderLayout.WEST);
		} else {
			this.getContentPane().remove(topRul);
			this.getContentPane().remove(leftRul);			
		}
		this.pack();
		this.repaint();
	}
	
	/** Set the underlying contact map data model */
	public void spawnNewViewWindow(Model mod) {
		String wintitle = "Contact Map of " + mod.getPDBCode() + " " + mod.getChainCode();
		View view = new View(mod, wintitle, Start.PYMOL_SERVER_URL);
		if(view == null) {
			System.err.println("Error: Couldn't initialize contact map window");
			//System.exit(-1);
		}
		System.out.println("Contact map " + mod.getPDBCode() + " " + mod.getChainCode() + " loaded.");

		// load structure in pymol
		view.pymolAdaptor = new PyMolAdaptor(this.pyMolServerUrl,	mod.getPDBCode(), mod.getChainCode(), mod.getTempPDBFileName());

		// if previous window was empty (not showing a contact map) dispose it
		if(this.mod == null) {
			this.setVisible(false);
			this.dispose();
		}
	}

	/** Returns the status variable currentAction which contains the currently selected actions */
	public int getCurrentAction(){
		return currentAction;
	}
	
//	/** Returns the current painting color set by the user */
//	public Color getCurrentPaintingColor() {
//		return currentPaintingColor;
//	}
	
	/** Returns the status variable doShowPdbSers which indicates whether PDB serials should be displayed */
	public boolean getShowPdbSers() {
		return doShowPdbSers;
	}

	/** Returns the status variable highlightComNbh which indicates whether cells 
	 * should be displayed in different colors based on common neighbourhoods sizes
	 */
	public boolean getHighlightComNbh() {
		return highlightComNbh;
	}

	public void setHighlightComNbh(boolean highlightComNbh) {
		this.highlightComNbh=highlightComNbh;
	}
	
	public HashMap<Contact,Integer> getComNbhSizes() {
		return comNbhSizes;
	}
}

