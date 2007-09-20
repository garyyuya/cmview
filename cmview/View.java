package cmview;
import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.awt.image.BufferedImage;
import javax.imageio.*;

import cmview.datasources.*;
import proteinstructure.*;

/**
 * Main GUI window and associated event handling.
 * Multiple instances of this will be shown in separate windows.
 * 
 * Initialized with mod=null, an empty window with the menu bars is shown.
 * Initialized with a valid model, the contact map is displayed in a ContactMapPane.
 * Multiple instances of this class are created to show multiple contact maps.
 */
public class View extends JFrame implements ActionListener {

	// constants
	static final long serialVersionUID = 1l;
	protected static final int SQUARE_SEL = 1;
	protected static final int FILL_SEL = 2;
	protected static final int NODE_NBH_SEL = 3;
	protected static final int SHOW_COMMON_NBH = 4;
	protected static final int DIAG_SEL = 5;
	protected static final int SEL_MODE_COLOR = 6;
	
	private static final boolean FIRST_MODEL = false;
	private static final boolean SECOND_MODEL = true;
	
	// menu item labels
	private static final String LABEL_FILE_INFO = "Info";
	private static final String LABEL_FILE_PRINT = "Print";	
	private static final String LABEL_FILE_QUIT = "Quit";	
	private static final String LABEL_DELETE_CONTACTS = "Delete selected contacts";
	private static final String LABEL_SHOW_TRIANGLES_3D = "Show Common Neighbour Triangles in 3D";
	private static final String LABEL_SHOW_COMMON_NBS_MODE = "Show Common Neighbours Mode";
	private static final String LABEL_SHOW_CONTACTS_3D = "Show Selected Contacts in 3D";
	private static final String LABEL_NODE_NBH_SELECTION_MODE = "Node Neighbourhood Selection Mode";
	private static final String LABEL_DIAGONAL_SELECTION_MODE = "Diagonal Selection Mode";
	private static final String LABEL_FILL_SELECTION_MODE = "Fill Selection Mode";
	private static final String LABEL_SQUARE_SELECTION_MODE = "Square Selection Mode";
	private static final String LABEL_SEL_MODE_COLOR = "Select by color mode";
	private static final String LABEL_COMPARE_CM = "Load second structure from"; 
	private static final String LABEL_SHOW_COMMON = "Show/hide common contacts";
	private static final String LABEL_SHOW_FIRST = "Show/hide contacts unique to first structure";
	private static final String LABEL_SHOW_SECOND = "Show/hide contacts unique to second structure";
	protected static final String LABEL_SHOW_PAIR_DIST_3D = "Show residue pair (%s,%s) as edge in 3D";	// used in ContactMapPane.showPopup
	
	// GUI components in the main frame
	JPanel statusPane; 			// panel holding the status bar (currently not used)
	JLabel statusBar; 			// TODO: Create a class StatusBar
	JToolBar toolBar;			// icon tool bar 
	JPanel cmp; 				// Main panel holding the Contact map pane
	JPanel topRul;				// Panel for top ruler	// TODO: Move this to ContactMapPane?
	JPanel leftRul;				// Panel for left ruler	// TODO: Move this to ContactMapPane?
	JPopupMenu popup; 			// right-click context menu
	JPanel tbPane;				// tool bar panel holding toolBar and cmp (necessary if toolbar is floatable)
	//JLayeredPane cmp2; 		// added for testing

	// Tool bar buttons
	JButton tbFileInfo, tbFilePrint, tbFileQuit, tbShowSel3D, tbShowComNbh3D,  tbDelete;  
	JToggleButton tbSquareSel, tbFillSel, tbDiagSel, tbNbhSel, tbShowComNbh, tbSelModeColor;
	JToggleButton tbViewPdbResSer, tbViewRuler, tbViewNbhSizeMap, tbViewDistanceMap, tbViewDensityMap, tbShowCommon, tbShowFirstStructure, tbShowSecondStructure;
	
	// Menu items
	JMenuItem sendM, squareM, fillM, comNeiM, triangleM, nodeNbhSelM, rangeM, delEdgesM, mmSelModeColor;
	JMenuItem sendP, squareP, fillP, comNeiP, triangleP, nodeNbhSelP, rangeP,  delEdgesP, popupSendEdge, pmSelModeColor;
	JMenuItem mmLoadGraph, mmLoadPdbase, mmLoadMsd, mmLoadCm, mmLoadPdb, mmLoadFtp;
	JMenuItem mmLoadGraph2, mmLoadPdbase2, mmLoadMsd2, mmLoadCm2, mmLoadPdb2, mmLoadFtp2;
	JMenuItem mmSaveGraphDb, mmSaveCmFile, mmSavePng;
	JMenuItem mmViewShowPdbResSers, mmViewHighlightComNbh, mmViewShowDensity, mmViewRulers, mmViewIconBar, mmViewShowDistMatrix;
	JMenuItem mmSelectAll, mmSelectByResNum, mmSelectHelixHelix, mmSelectBetaBeta, mmSelectInterSsContacts, mmSelectIntraSsContacts;
	JMenuItem mmColorReset, mmColorPaint, mmColorChoose;
	JMenuItem mmSelCommonContactsInComparedMode,  mmSelFirstStrucInComparedMode,  mmSelSecondStrucInComparedMode;
	JMenuItem mmToggleDiffDistMap;
	JMenuItem mmInfo, mmPrint, mmQuit, mmHelpAbout, mmHelpHelp, mmHelpWriteConfig;
	
	// Data and status variables
	private Model mod;
	public ContactMapPane cmPane;
	public ResidueRuler topRuler;
	public ResidueRuler leftRuler;
	private int pymolSelSerial;		 	// for incremental numbering // TODO: Move this to PymolAdaptor
	private int pymolNbhSerial;			// for incremental numbering // TODO: Move this to PymolAdaptor
	
	private boolean common;
	private boolean firstS;
	private boolean secondS;

	// current gui state
	private int currentSelectionMode;	// current selection mode (see constants above), modify using setSelectionMode
	private boolean showPdbSers;		// whether showing pdb serials is switched on
	private boolean showRulers;			// whether showing residue rulers is switched on
	private boolean showIconBar;		// whether showing the icon bar is switched on
	private boolean showNbhSizeMap;		// whether showing the common neighbourhood size map is switched on 
	private boolean showDensityMap;		// whether showing the density map is switched on
	private boolean showDistanceMap;	// whether showing the distance map is switched on
	private boolean compareStatus;		// tells ContactMapPane to draw compared contact map if 2. structure is loaded
	private Color currentPaintingColor;	// current color for coloring contacts selected by the user
	private boolean selCommonContactsInComparedMode;	// when true, selection on compared contact map in both structures possible
	private boolean selFirstStrucInComparedMode; // when true selection on compared contact map in first structure possible
	private boolean selSecondStrucInComparedMode; // when true selection on compared contact map in second structure possible
	private boolean showDiffDistMap; 	// whether showing the difference distance map is switched on
	private boolean comparisonMode;		// whether comparison functions are enabled
	
	// global icons TODO: replace these by tickboxes?
	ImageIcon icon_selected = new ImageIcon(this.getClass().getResource("/icons/tick.png"));
	ImageIcon icon_deselected = new ImageIcon(this.getClass().getResource("/icons/bullet_blue.png"));

	/*----------------------------- constructors ----------------------------*/
	
	/** Create a new View object */
	public View(Model mod, String title) {
		super(title);
		this.mod = mod;
		if(mod == null) {
			this.setPreferredSize(new Dimension(Start.INITIAL_SCREEN_SIZE,Start.INITIAL_SCREEN_SIZE));
		}
		this.currentSelectionMode = SQUARE_SEL;
		this.pymolSelSerial = 1;
		this.pymolNbhSerial = 1;
		this.showPdbSers = false;
		this.showNbhSizeMap = false;
		this.showRulers=Start.SHOW_RULERS_ON_STARTUP;
		this.showIconBar=Start.SHOW_ICON_BAR;
		this.showDensityMap=false;
		this.showDistanceMap=false;
		this.currentPaintingColor = Color.blue;
		this.selCommonContactsInComparedMode= true;
		this.selFirstStrucInComparedMode= true;
		this.selSecondStrucInComparedMode= true;
		this.showDiffDistMap = false;
		this.comparisonMode = false;
		
		this.initGUI(); // build gui tree and show window
		this.toFront(); // bring window to front
		this.compareStatus = false;	
	}
	
	/*---------------------------- private methods --------------------------*/
	
	/**
	 * Sets the current selection mode. This sets the internal state variable and changes some gui components.
	 * Use getSelectionMode to retrieve the current state.
	 */
	private void setSelectionMode(int mode) {
		switch(mode) {
		case(SQUARE_SEL): tbSquareSel.setSelected(true); break;
		case(FILL_SEL): tbFillSel.setSelected(true); break;
		case(NODE_NBH_SEL): tbNbhSel.setSelected(true); break;
		case(SHOW_COMMON_NBH): tbShowComNbh.setSelected(true); break;
		case(DIAG_SEL): tbDiagSel.setSelected(true); break;
		case(SEL_MODE_COLOR) : tbSelModeColor.setSelected(true); break;
		default: System.err.println("Error in setSelectionMode. Unknown selection mode " + mode); return;
		}
		this.currentSelectionMode = mode;
	}
	
	/**
	 * Sets up and returns a new menu item with the given icon and label, adds it to the given JMenu and
	 * registers 'this' as the action listener.
	 */
	private JMenuItem makeMenuItem(String label, Icon icon, JMenu menu) {
		JMenuItem newItem = new JMenuItem(label, icon);
		newItem.addActionListener(this);
		menu.add(newItem);
		return newItem;
	}
	
	/**
	 * Sets up and returns a new popup menu item with the given icon and label, adds it to the given JPopupMenu and
	 * registers 'this' as the action listener.
	 */	
	private JMenuItem makePopupMenuItem(String label, Icon icon, JPopupMenu menu) {
		JMenuItem newItem = new JMenuItem(label, icon);
		newItem.addActionListener(this);
		menu.add(newItem);
		return newItem;		
	}

	/**
	 * Sets up and returns a new tool bar button
	 */
	private JButton makeToolBarButton(ImageIcon icon, String toolTipText) {
		JButton newButton = new JButton(icon);
		newButton.setFocusPainted(false);
		newButton.setToolTipText(toolTipText);
		newButton.addActionListener(this);
		toolBar.add(newButton);
		return newButton;
	}
	
	/**
	 * Sets up and returns a new tool bar toggle button
	 */
	private JToggleButton makeToolBarToggleButton(ImageIcon icon, String toolTipText, boolean selected, boolean enabled, boolean visible) {
		JToggleButton newButton = new JToggleButton(icon, selected);
		newButton.setFocusPainted(false);
		newButton.setToolTipText(toolTipText);
		newButton.addActionListener(this);
		newButton.setVisible(visible);
		newButton.setEnabled(enabled);
		newButton.setSelected(selected);		
		toolBar.add(newButton);
		return newButton;
	}
	
	/**
	 * Creates an icon which changes the color with the variable 'currentPaintingColor'.
	 * @return The 'magic' icon
	 */
	private Icon getCurrentColorIcon() {
		return new Icon() {
			public void paintIcon(Component c, Graphics g, int x, int y) {
				Color oldColor = c.getForeground();
				g.setColor(currentPaintingColor);
				g.translate(x, y);
				g.fillRect(2,2,12,12);
				g.translate(-x, -y);  // Restore Graphics object
				g.setColor(oldColor);
			}			
			public int getIconHeight() {
				return 16;
			}
			public int getIconWidth() {
				return 16;
			}
		};
	}
	
	/**
	 * Creates a black-square icon.
	 * @return The black icon
	 */
	private Icon getBlackSquareIcon() {
		return new Icon() {
			public void paintIcon(Component c, Graphics g, int x, int y) {
				Color oldColor = c.getForeground();
				g.setColor(Color.black);
				g.translate(x, y);
				g.fillRect(2,2,12,12);
				g.translate(-x, -y);  // Restore Graphics object
				g.setColor(oldColor);
			}	
			public int getIconHeight() {
				return 16;
			}
			public int getIconWidth() {
				return 16;
			}
		};
	}
	
	/** Initialize and show the main GUI window */
	private void initGUI(){

		// Setting the main layout 
		setLayout(new BorderLayout());
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setLocation(20,20);

		// Creating the Panels
		tbPane = new JPanel(new BorderLayout());	// toolbar pane, holding only toolbar and cmp (all the rest)
		cmp = new JPanel(new BorderLayout()); 		// pane holding the cmPane and (optionally) the ruler panes
		topRul = new JPanel(new BorderLayout()); 	// pane holding the top ruler
		leftRul = new JPanel(new BorderLayout()); 	// pane holding the left ruler
		statusPane = new JPanel(); 					// pane holding the status bar, TODO: Create a class StatusBar
		statusBar = new JLabel(" ");				// a primitive status bar for testing
		//cmp2 = new JLayeredPane(); 				// for testing with layered panes
								
		// Icons
		ImageIcon icon_square_sel_mode = new ImageIcon(this.getClass().getResource("/icons/shape_square.png"));
		ImageIcon icon_fill_sel_mode = new ImageIcon(this.getClass().getResource("/icons/paintcan.png"));
		ImageIcon icon_diag_sel_mode = new ImageIcon(this.getClass().getResource("/icons/diagonals.png"));
		ImageIcon icon_nbh_sel_mode = new ImageIcon(this.getClass().getResource("/icons/group.png"));
		ImageIcon icon_show_sel_cont_3d = new ImageIcon(this.getClass().getResource("/icons/shape_square_go.png"));
		ImageIcon icon_show_com_nbs_mode = new ImageIcon(this.getClass().getResource("/icons/shape_flip_horizontal.png"));
		ImageIcon icon_show_triangles_3d = new ImageIcon(this.getClass().getResource("/icons/shape_rotate_clockwise.png"));
		ImageIcon icon_del_contacts = new ImageIcon(this.getClass().getResource("/icons/cross.png"));	
		ImageIcon icon_show_pair_dist_3d = new ImageIcon(this.getClass().getResource("/icons/user_go.png"));
		ImageIcon icon_colorwheel = new ImageIcon(this.getClass().getResource("/icons/color_wheel.png"));
		ImageIcon icon_file_info = new ImageIcon(this.getClass().getResource("/icons/information.png"));
		ImageIcon icon_file_print = new ImageIcon(this.getClass().getResource("/icons/printer.png"));		
		ImageIcon icon_file_quit = new ImageIcon(this.getClass().getResource("/icons/door_open.png"));				
		ImageIcon icon_show_common = new ImageIcon(this.getClass().getResource("/icons/page_copy.png"));
		ImageIcon icon_show_first = new ImageIcon(this.getClass().getResource("/icons/page_delete.png"));
		ImageIcon icon_show_second = new ImageIcon(this.getClass().getResource("/icons/page_add.png"));
		ImageIcon icon_sel_mode_color = new ImageIcon(this.getClass().getResource("/icons/color_swatch.png"));
		Icon icon_color = getCurrentColorIcon();	// magic icon with current painting color
		Icon icon_black = getBlackSquareIcon();		// black square icon
		
		// Tool bar
		toolBar = new JToolBar();
		toolBar.setVisible(Start.SHOW_ICON_BAR);
		
		tbFileInfo = makeToolBarButton(icon_file_info, LABEL_FILE_INFO);
		tbFilePrint = makeToolBarButton(icon_file_print, LABEL_FILE_PRINT);
		tbFileQuit = makeToolBarButton(icon_file_quit, LABEL_FILE_QUIT);
		toolBar.addSeparator();
		tbSquareSel = makeToolBarToggleButton(icon_square_sel_mode, LABEL_SQUARE_SELECTION_MODE, true, true, true);
		tbFillSel = makeToolBarToggleButton(icon_fill_sel_mode, LABEL_FILL_SELECTION_MODE, false, true, true);
		tbDiagSel = makeToolBarToggleButton(icon_diag_sel_mode, LABEL_DIAGONAL_SELECTION_MODE, false, true, true);
		tbNbhSel = makeToolBarToggleButton(icon_nbh_sel_mode, LABEL_NODE_NBH_SELECTION_MODE, false, true, true);
		tbShowComNbh = makeToolBarToggleButton(icon_show_com_nbs_mode, LABEL_SHOW_COMMON_NBS_MODE, false, true, true);
		tbSelModeColor = makeToolBarToggleButton(icon_sel_mode_color, LABEL_SEL_MODE_COLOR, false, true, true);
		toolBar.addSeparator();		
		tbShowSel3D = makeToolBarButton(icon_show_sel_cont_3d, LABEL_SHOW_CONTACTS_3D);
		tbShowComNbh3D = makeToolBarButton(icon_show_triangles_3d, LABEL_SHOW_TRIANGLES_3D);
		tbDelete = makeToolBarButton(icon_del_contacts, LABEL_DELETE_CONTACTS);
		toolBar.addSeparator(new Dimension(100, 10));
		tbShowCommon = makeToolBarToggleButton(icon_show_common, LABEL_SHOW_COMMON, selCommonContactsInComparedMode, true, false);
		tbShowFirstStructure = makeToolBarToggleButton(icon_show_first, LABEL_SHOW_FIRST, selFirstStrucInComparedMode, true, false);
		tbShowSecondStructure = makeToolBarToggleButton(icon_show_second, LABEL_SHOW_SECOND, selSecondStrucInComparedMode, true, false);
		
		// Toggle buttons in view menu (not being used yet)
		tbViewPdbResSer = new JToggleButton();
		tbViewRuler = new JToggleButton();
		tbViewNbhSizeMap = new JToggleButton();
		tbViewDistanceMap = new JToggleButton();
		tbViewDensityMap = new JToggleButton();
		toolBar.setFloatable(Start.ICON_BAR_FLOATABLE);
		
		// ButtonGroup for selection modes (so upon selecting one, others are deselected automatically)
		ButtonGroup selectionModeButtons = new ButtonGroup();
		selectionModeButtons.add(tbSquareSel);
		selectionModeButtons.add(tbFillSel);
		selectionModeButtons.add(tbDiagSel);
		selectionModeButtons.add(tbSelModeColor);
		selectionModeButtons.add(tbNbhSel);
		selectionModeButtons.add(tbShowComNbh);
		
		// Popup menu
		JPopupMenu.setDefaultLightWeightPopupEnabled(false);
		popup = new JPopupMenu();
		
		squareP = makePopupMenuItem(LABEL_SQUARE_SELECTION_MODE, icon_square_sel_mode, popup);
		fillP = makePopupMenuItem(LABEL_FILL_SELECTION_MODE, icon_fill_sel_mode, popup);
		rangeP = makePopupMenuItem(LABEL_DIAGONAL_SELECTION_MODE, icon_diag_sel_mode, popup);
		nodeNbhSelP = makePopupMenuItem(LABEL_NODE_NBH_SELECTION_MODE, icon_nbh_sel_mode, popup);
		comNeiP = makePopupMenuItem(LABEL_SHOW_COMMON_NBS_MODE, icon_show_com_nbs_mode, popup);
		pmSelModeColor = makePopupMenuItem(LABEL_SEL_MODE_COLOR, icon_sel_mode_color, popup);
		if (Start.USE_PYMOL) {
			popup.addSeparator();		
			sendP = makePopupMenuItem(LABEL_SHOW_CONTACTS_3D, icon_show_sel_cont_3d, popup);
			popupSendEdge = makePopupMenuItem(LABEL_SHOW_PAIR_DIST_3D, icon_show_pair_dist_3d, popup);
			triangleP = makePopupMenuItem(LABEL_SHOW_TRIANGLES_3D, icon_show_triangles_3d, popup);
		}
		popup.addSeparator();		
		delEdgesP = makePopupMenuItem(LABEL_DELETE_CONTACTS, icon_del_contacts, popup);

		// Main menu
		JMenuBar menuBar;
		JMenu menu, submenu;
		menuBar = new JMenuBar();
		this.setJMenuBar(menuBar);

		// File menu
		menu = new JMenu("File");
		menu.setMnemonic(KeyEvent.VK_F);
		mmInfo = makeMenuItem(LABEL_FILE_INFO, null, menu);		
		// Load
		submenu = new JMenu("Load from");
		if(Start.USE_DATABASE) {
			mmLoadGraph = makeMenuItem("Graph database",null,submenu);
			mmLoadPdbase = makeMenuItem("Pdbase",null,submenu);
			mmLoadMsd = makeMenuItem("MSD",null, submenu);
		}		
		mmLoadFtp = makeMenuItem("Online PDB", null, submenu);
		mmLoadPdb = makeMenuItem("PDB file", null, submenu);
		mmLoadCm = makeMenuItem("Contact map file", null, submenu);		
		menu.add(submenu);
		// Save
		submenu = new JMenu("Save to");
		mmSaveCmFile = makeMenuItem("Contact map file", null, submenu);
		mmSavePng = makeMenuItem("PNG file", null, submenu);
		if(Start.USE_DATABASE) {
			mmSaveGraphDb = makeMenuItem("Graph database", null, submenu);
		}
		menu.add(submenu);		
		// Print, Quit
		mmPrint = makeMenuItem(LABEL_FILE_PRINT, null, menu);
		mmQuit = makeMenuItem(LABEL_FILE_QUIT, null, menu);
		menuBar.add(menu);

		// View menu
		menu = new JMenu("View");
		menu.setMnemonic(KeyEvent.VK_V);		
		mmViewShowPdbResSers = makeMenuItem("Toggle show PDB residue numbers", icon_deselected, menu);
		mmViewRulers = makeMenuItem("Toggle rulers", icon_deselected, menu);
		mmViewIconBar = makeMenuItem("Toggle icon bar", icon_deselected, menu);	// doesn't work properly if icon bar is floatable
		menu.addSeparator();		
		mmViewHighlightComNbh = makeMenuItem("Toggle highlight of cells by common neighbourhood size", icon_deselected, menu);
		mmViewShowDensity = makeMenuItem("Toggle show contact density", icon_deselected, menu);
		mmViewShowDistMatrix = makeMenuItem("Toggle show distance matrix", icon_deselected, menu);
		menuBar.add(menu);

		// Select menu
		menu = new JMenu("Select");
		menu.setMnemonic(KeyEvent.VK_S);
		mmSelectAll = makeMenuItem("All contacts", null, menu);
		mmSelectByResNum = makeMenuItem("By residue number", null, menu);
		menu.addSeparator();
		mmSelectHelixHelix = makeMenuItem("Helix-Helix contacts", null, menu);
		mmSelectBetaBeta = makeMenuItem("Strand-Strand contacts", null, menu);
		mmSelectInterSsContacts = makeMenuItem("Contacts between SS elements", null, menu);
		mmSelectIntraSsContacts = makeMenuItem("Contacts within SS elements", null, menu);		
		menuBar.add(menu);
		
		// Color menu
		menu = new JMenu("Color");
		menu.setMnemonic(KeyEvent.VK_C);
		mmColorChoose = makeMenuItem("Choose painting color", icon_colorwheel, menu);
		mmColorPaint = makeMenuItem("Color selected contacts", icon_color, menu);
		mmColorReset= makeMenuItem("Reset contact colors to black", icon_black, menu);
		menuBar.add(menu);
		
		// Action menu, TODO: split into 'Mode' and 'Action'
		menu = new JMenu("Action");
		menu.setMnemonic(KeyEvent.VK_A);
		squareM = makeMenuItem(LABEL_SQUARE_SELECTION_MODE, icon_square_sel_mode, menu);
		fillM = makeMenuItem(LABEL_FILL_SELECTION_MODE, icon_fill_sel_mode, menu);
		rangeM = makeMenuItem(LABEL_DIAGONAL_SELECTION_MODE,icon_diag_sel_mode, menu);
		nodeNbhSelM = makeMenuItem(LABEL_NODE_NBH_SELECTION_MODE, icon_nbh_sel_mode, menu);
		comNeiM = makeMenuItem(LABEL_SHOW_COMMON_NBS_MODE, icon_show_com_nbs_mode, menu);
		mmSelModeColor = makeMenuItem(LABEL_SEL_MODE_COLOR, icon_sel_mode_color, menu);
		if (Start.USE_PYMOL) {
			menu.addSeparator();
			sendM = makeMenuItem(LABEL_SHOW_CONTACTS_3D, icon_show_sel_cont_3d, menu);
			triangleM = makeMenuItem(LABEL_SHOW_TRIANGLES_3D, icon_show_triangles_3d, menu);
		}
		menu.addSeparator();		
		delEdgesM = makeMenuItem(LABEL_DELETE_CONTACTS, icon_del_contacts, menu);
		menuBar.add(menu);

		// Comparison Menu
		menu = new JMenu("Compare");
		menu.setMnemonic(KeyEvent.VK_P);
		// Load
		submenu = new JMenu(LABEL_COMPARE_CM);
		menu.add(submenu);		
		if(Start.USE_DATABASE) {
			mmLoadGraph2 = makeMenuItem("Graph database",null,submenu);
			mmLoadPdbase2 = makeMenuItem("Pdbase",null,submenu);
			mmLoadMsd2 = makeMenuItem("MSD",null, submenu);
		}		
		mmLoadFtp2 = makeMenuItem("Online PDB", null, submenu);
		mmLoadPdb2 = makeMenuItem("PDB file", null, submenu);
		mmLoadCm2 = makeMenuItem("Contact map file", null, submenu);		
		menu.addSeparator();
		mmSelCommonContactsInComparedMode = makeMenuItem("Toggle show common contacts", icon_selected, menu);
		mmSelFirstStrucInComparedMode = makeMenuItem("Toggle show contacts unique in first structure", icon_selected, menu);
		mmSelSecondStrucInComparedMode = makeMenuItem("Toggle show contacts unique in second structure ", icon_selected, menu);		
		menu.addSeparator();
		mmToggleDiffDistMap = makeMenuItem("Toggle show difference map", icon_deselected, menu);
		menuBar.add(menu);
		
		// Help menu
		menu = new JMenu("Help");
		menu.setMnemonic(KeyEvent.VK_H);	
		mmHelpHelp = makeMenuItem("Help", null, menu);
		mmHelpWriteConfig = makeMenuItem("Write example configuration file", null, menu);
		mmHelpAbout = makeMenuItem("About", null, menu);
		menuBar.add(menu);
		
		// Status bar
//		if (mod!= null){
//			int[] statusRuler = new int[3];
//			statusRuler = cmPane.getRulerCoordinates();
//			String r1 = "" + statusRuler[0] + "";
//			String r2 = "" + statusRuler[1] + "";
//			String r3 = "" + statusRuler[2] + "";
//
//
//			statusBar.setText(r1 +" , " + r2 + " , " +r3);
//			statusPane.setLayout(new BorderLayout());
//			statusPane.add(statusBar, BorderLayout.CENTER);
//			statusPane.setBorder(new EtchedBorder(EtchedBorder.LOWERED));
//		}
		
		// Creating contact map pane if model loaded
		if(mod != null) {
			// testing: use two separate layers for contacts and interaction
//			TestCmBackgroundPane cmPaneBg = new TestCmBackgroundPane(mod, this);
//			TestCmInteractionPane cmPaneInt = new TestCmInteractionPane(mod, this, cmPaneBg);
//			cmp.setLayout(new OverlayLayout(cmp));
//			cmp.add(cmPaneInt);
//			cmp.add(cmPaneBg);

			// alternative: use JLayeredPane
//			TestCmInteractionPane cmPaneInt = new TestCmInteractionPane(mod, this);
//			TestCmBackgroundPane cmPaneBg = new TestCmBackgroundPane(mod, this);
//			cmp2.setLayout(new OverlayLayout(cmp2));
//			cmp2.add(cmPaneInt, new Integer(2));
//			cmp2.add(cmPaneBg, new Integer(1));
			 //add cmp2 to contentPane
			
			cmPane = new ContactMapPane(mod, this);
			cmp.add(cmPane, BorderLayout.CENTER);
			topRuler = new ResidueRuler(cmPane,mod,this,ResidueRuler.TOP);
			leftRuler = new ResidueRuler(cmPane,mod,this,ResidueRuler.LEFT);
			topRul.add(topRuler);
			leftRul.add(leftRuler);
		}
		
		// Add everything to the content pane		
		this.tbPane.add(toolBar, BorderLayout.NORTH);			// tbPane is necessary if toolBar is floatable
		this.tbPane.add(cmp,BorderLayout.CENTER);				// otherwise can add these to contentPane directly
		this.getContentPane().add(tbPane, BorderLayout.CENTER); // and get rid of this line
		if(showRulers) {
			cmp.add(topRul, BorderLayout.NORTH);
			cmp.add(leftRul, BorderLayout.WEST);
		}
		//this.getContentPane().add(statusPane,BorderLayout.SOUTH);

		// Show GUI
		pack();
		setVisible(true);
	}

	/*------------------------------ event handling -------------------------*/

	/**
	 * Handling action events for all menu items.
	 */
	public void actionPerformed (ActionEvent e) {

		/* ---------- File Menu ---------- */
		
		// Load
		if(e.getSource() == mmLoadGraph) {
			handleLoadFromGraphDb(FIRST_MODEL);
		}
		if(e.getSource() == mmLoadPdbase) {
			handleLoadFromPdbase(FIRST_MODEL);

		}		  
		if(e.getSource() == mmLoadMsd) {
			handleLoadFromMsd(FIRST_MODEL);
		}			  
		if(e.getSource() == mmLoadPdb) {
			handleLoadFromPdbFile(FIRST_MODEL);
		}
		if(e.getSource() == mmLoadFtp) {
			handleLoadFromFtp(FIRST_MODEL);
		}
		if(e.getSource() == mmLoadCm) {
			handleLoadFromCmFile(FIRST_MODEL);
		}
		
		// Save
		if(e.getSource() == mmSaveGraphDb) {
			handleSaveToGraphDb();
		}		  
		if(e.getSource() == mmSaveCmFile) {
			handleSaveToCmFile();
		}
		if(e.getSource() == mmSavePng) {
			handleSaveToPng();
		}	
		
		// Info, Print, Quit
		if(e.getSource() == mmInfo || e.getSource() == tbFileInfo) {
			handleInfo();
		}		  		  
		if(e.getSource() == mmPrint || e.getSource() == tbFilePrint) {
			handlePrint();
		}		  
		if(e.getSource() == mmQuit || e.getSource() == tbFileQuit) {
			handleQuit();
		}
		
		/* ---------- View Menu ---------- */		
		
		if(e.getSource() == mmViewShowPdbResSers) {
			handleShowPdbResSers();
		}		  
		if(e.getSource() == mmViewRulers) {
			handleShowRulers();
		}	
		if(e.getSource() == mmViewIconBar) {
			handleShowIconBar();
		}	
		if(e.getSource() == mmViewHighlightComNbh) {
			handleShowNbhSizeMap();
		}
		if(e.getSource() == mmViewShowDensity) {
			handleShowDensityMap();
		}
		if(e.getSource() == mmViewShowDistMatrix) {
			handleShowDistanceMap();
		}
		
		/* ---------- Select menu ---------- */
		
		if(e.getSource() == mmSelectAll  ) {
			handleSelectAll();
		}
		if(e.getSource() == mmSelectByResNum  ) {
			handleSelectByResNum();
		}
		if(e.getSource() == mmSelectHelixHelix  ) {
			handleSelectHelixHelix();
		}
		if(e.getSource() == mmSelectBetaBeta  ) {
			handleSelectBetaBeta();
		}
		if(e.getSource() == mmSelectInterSsContacts  ) {
			handleSelectInterSsContacts();
		}
		if(e.getSource() == mmSelectIntraSsContacts  ) {
			handleSelectIntraSsContacts();
		}
		
		/* ---------- Color menu ---------- */
		
		if(e.getSource() == mmColorReset) {
			handleColorReset();
		}
		if(e.getSource() == mmColorPaint) {
			handleColorPaint();
		}
		if(e.getSource() == mmColorChoose) {
			handleColorSelect();
		}

		/* ---------- Action menu ---------- */

		// Selection modes
		
		// square button clicked
		if (e.getSource() == squareM || e.getSource() == squareP || e.getSource() == tbSquareSel ) {
			setSelectionMode(SQUARE_SEL);
		}
		// fill button clicked
		if (e.getSource() == fillM || e.getSource() == fillP || e.getSource() == tbFillSel) {
			setSelectionMode(FILL_SEL);
		}			
		// diagonal selection clicked
		if (e.getSource() == rangeM || e.getSource() == rangeP || e.getSource() == tbDiagSel) {
			setSelectionMode(DIAG_SEL);
		}
		// node neihbourhood selection button clicked 
		if (e.getSource() == nodeNbhSelM || e.getSource() == nodeNbhSelP || e.getSource() == tbNbhSel) {
			setSelectionMode(NODE_NBH_SEL);
		}		
		// showing com. Nei. button clicked
		if (e.getSource() == comNeiM || e.getSource() == comNeiP || e.getSource() == tbShowComNbh) {
			setSelectionMode(SHOW_COMMON_NBH);
		}
		// color selection mode button clicked
		if (e.getSource() == mmSelModeColor || e.getSource() == pmSelModeColor || e.getSource() == tbSelModeColor) {
			setSelectionMode(SEL_MODE_COLOR);
		}		
		
		// Actions
		
		// send selection button clicked
		if (e.getSource() == sendM || e.getSource() == sendP || e.getSource() == tbShowSel3D) {	
			handleShowSelContacts3D();
		}
		// send com.Nei. button clicked
		if(e.getSource()== triangleM || e.getSource()== triangleP || e.getSource() == tbShowComNbh3D) {
			handleShowTriangles3D();
		}		
		// delete selected edges button clicked
		if (e.getSource() == delEdgesM || e.getSource() == delEdgesP || e.getSource() == tbDelete) {
			handleDeleteSelContacts();
		}
		// send current edge (only available in popup menu)
		if(e.getSource() == popupSendEdge) {
			handleShowDistance3D();
		}
		
		/* ---------- Comparison Menu ---------- */
		
		
		/** for the contact map comparison load menu */
		
		
		if(e.getSource() == mmLoadGraph2) {
			handleLoadFromGraphDb(SECOND_MODEL);
		}
		if(e.getSource() == mmLoadPdbase2) {
			handleLoadFromPdbase(SECOND_MODEL);
		}		  
		if(e.getSource() == mmLoadMsd2) {
			handleLoadFromMsd(SECOND_MODEL);
		}			  
		if(e.getSource() == mmLoadPdb2) {
			handleLoadFromPdbFile(SECOND_MODEL);
		}
		if(e.getSource() == mmLoadFtp2) {
			handleLoadFromFtp(SECOND_MODEL);
		}
		if(e.getSource() == mmLoadCm2) {
			handleLoadFromCmFile(SECOND_MODEL);
		}
		
		if(e.getSource() == mmSelCommonContactsInComparedMode || e.getSource() == tbShowCommon) {
			handleSelContactsInComparedMode();
		}
		
		if(e.getSource() == mmSelFirstStrucInComparedMode || e.getSource() == tbShowFirstStructure) {
			handleSelFirstStrucInComparedMode();
		}
		
		if(e.getSource() == mmSelSecondStrucInComparedMode || e.getSource() == tbShowSecondStructure) {
			handleSelSecondStrucInComparedMode();
		}
		
		if(e.getSource() == mmToggleDiffDistMap) {
			handleToggleDiffDistMap();
		}
		
		/* ---------- Help Menu ---------- */
		
		if(e.getSource() == mmHelpAbout) {
			handleHelpAbout();
		}
		if(e.getSource() == mmHelpHelp) {
			handleHelpHelp();
		}
		if(e.getSource() == mmHelpWriteConfig) {
			handleHelpWriteConfig();
		}
		
	}
	

	/* -------------------- File Menu -------------------- */
	
	private void handleLoadFromGraphDb(boolean secondModel) {
		if(!Start.isDatabaseConnectionAvailable()) {
			showNoDatabaseConnectionWarning();
		} else {
			if (secondModel == SECOND_MODEL && mod == null){
				this.showNoContactMapWarning();
			} else{
				LoadDialog dialog = new LoadDialog(this, "Load from graph database", new LoadAction(secondModel) {
					public void doit(Object o, String f, String ac, String cc, String ct, double dist, int minss, int maxss, String db, int gid) {
						View view = (View) o;
						view.doLoadFromGraphDb(db, gid, secondModel);
					}
				}, null, null, null, null, null, null, null, Start.DEFAULT_GRAPH_DB, "");
				dialog.showIt();
			}
		}
	}

	public void doLoadFromGraphDb(String db, int gid, boolean secondModel) {
		System.out.println("Loading from graph database");
		System.out.println("Database:\t" + db);
		System.out.println("Graph Id:\t" + gid);
		try {
			Model mod = new GraphDbModel(gid, db);
			if(secondModel) {
				if(cmPane.addSecondModel(mod) == false) {
					showDifferentSizeContactMapsError();
				}else {
					compareStatus = true;
					cmPane.toggleCompareMode(compareStatus);
					Start.getPyMolAdaptor().loadStructure(mod.getTempPdbFileName(), mod.getPDBCode(), mod.getChainCode(), true);
					Start.getPyMolAdaptor().alignStructure(cmPane.getFirstModel().getPDBCode(), cmPane.getFirstModel().getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
					
				}
			} else {
				this.spawnNewViewWindow(mod);
			}
		} catch(ModelConstructionError e) {
			showLoadError(e.getMessage());
		}
	}

	private void handleLoadFromPdbase(boolean secondModel) {
		if(!Start.isDatabaseConnectionAvailable()) {
			showNoDatabaseConnectionWarning();
		} else {
			if (secondModel == SECOND_MODEL && mod == null){
				this.showNoContactMapWarning();
			} else{
				LoadDialog dialog = new LoadDialog(this, "Load from Pdbase", new LoadAction(secondModel) {
					public void doit(Object o, String f, String ac, String cc, String ct, double dist, int minss, int maxss, String db, int gid) {
						View view = (View) o;
						view.doLoadFromPdbase(ac, cc, ct, dist, minss, maxss, db, secondModel);
					}
				}, null, "", "", Start.DEFAULT_CONTACT_TYPE, String.valueOf(Start.DEFAULT_DISTANCE_CUTOFF), "", "", Start.DEFAULT_PDB_DB, null);
				dialog.showIt();
			}
		}
	}

	public void doLoadFromPdbase(String ac, String cc, String ct, double dist, int minss, int maxss, String db, boolean secondModel) {
		System.out.println("Loading from Pdbase");
		System.out.println("PDB code:\t" + ac);
		System.out.println("Chain code:\t" + cc);
		System.out.println("Contact type:\t" + ct);
		System.out.println("Dist. cutoff:\t" + dist);	
		System.out.println("Min. Seq. Sep.:\t" + (minss==-1?"none":minss));
		System.out.println("Max. Seq. Sep.:\t" + (maxss==-1?"none":maxss));
		System.out.println("Database:\t" + db);	
		try{
			Model mod = new PdbaseModel(ac, cc, ct, dist, minss, maxss, db);
			if(secondModel) {
				if(cmPane.addSecondModel(mod) == false) {
					showDifferentSizeContactMapsError();
				} else {
					compareStatus = true;
					cmPane.toggleCompareMode(compareStatus);
					Start.getPyMolAdaptor().loadStructure(mod.getTempPdbFileName(), mod.getPDBCode(), mod.getChainCode(), true);
					Start.getPyMolAdaptor().alignStructure(cmPane.getFirstModel().getPDBCode(), cmPane.getFirstModel().getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
					
				}
			 
			} else {
				this.spawnNewViewWindow(mod);
			}
		} catch(ModelConstructionError e) {
			showLoadError(e.getMessage());
		}
	}

	private void handleLoadFromMsd(boolean secondModel) {
		if(!Start.isDatabaseConnectionAvailable()) {
			showNoDatabaseConnectionWarning();
		} else {
			if (secondModel == SECOND_MODEL && mod == null){
				this.showNoContactMapWarning();
			} else{
				LoadDialog dialog = new LoadDialog(this, "Load from MSD", new LoadAction(secondModel) {
					public void doit(Object o, String f, String ac, String cc, String ct, double dist, int minss, int maxss, String db, int gid) {
						View view = (View) o;
						view.doLoadFromMsd(ac, cc, ct, dist, minss, maxss, db, secondModel);
					}
				}, null, "", "", Start.DEFAULT_CONTACT_TYPE, String.valueOf(Start.DEFAULT_DISTANCE_CUTOFF), "", "", Start.DEFAULT_MSDSD_DB, null);
				dialog.showIt();
			}
		}
	}

	public void doLoadFromMsd(String ac, String cc, String ct, double dist, int minss, int maxss, String db, boolean secondModel) {
		System.out.println("Loading from MSD");
		System.out.println("PDB code:\t" + ac);
		System.out.println("Chain code:\t" + cc);
		System.out.println("Contact type:\t" + ct);
		System.out.println("Dist. cutoff:\t" + dist);	
		System.out.println("Min. Seq. Sep.:\t" + (minss==-1?"none":minss));
		System.out.println("Max. Seq. Sep.:\t" + (maxss==-1?"none":maxss));
		System.out.println("Database:\t" + db);	
		try {
			Model mod = new MsdsdModel(ac, cc, ct, dist, minss, maxss, db);
			if(secondModel) {
				if(cmPane.addSecondModel(mod) == false) {
					showDifferentSizeContactMapsError();
				}else {
					compareStatus = true;
					cmPane.toggleCompareMode(compareStatus);
					Start.getPyMolAdaptor().loadStructure(mod.getTempPdbFileName(), mod.getPDBCode(), mod.getChainCode(), true);
					Start.getPyMolAdaptor().alignStructure(cmPane.getFirstModel().getPDBCode(), cmPane.getFirstModel().getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
//					disable all old objects and selections
				}
			} else {
				this.spawnNewViewWindow(mod);
			}
		} catch(ModelConstructionError e) {
			showLoadError(e.getMessage());
		}
	}	  

	private void handleLoadFromPdbFile(boolean secondModel) {
		
		if (secondModel == SECOND_MODEL && mod == null){
			this.showNoContactMapWarning();
		} else{
			LoadDialog dialog = new LoadDialog(this, "Load from Pdb file", new LoadAction(secondModel) {
				public void doit(Object o, String f, String ac, String cc, String ct, double dist, int minss, int maxss, String db, int gid) {
					View view = (View) o;
					view.doLoadFromPdbFile(f, cc, ct, dist, minss, maxss, secondModel);
				}
			}, "", null, "", Start.DEFAULT_CONTACT_TYPE, String.valueOf(Start.DEFAULT_DISTANCE_CUTOFF), "", "", null, null);
			dialog.showIt();
		}
	}

	public void doLoadFromPdbFile(String f, String cc, String ct, double dist, int minss, int maxss, boolean secondModel) {
		System.out.println("Loading from Pdb file");
		System.out.println("Filename:\t" + f);
		System.out.println("Chain code:\t" + cc);
		System.out.println("Contact type:\t" + ct);
		System.out.println("Dist. cutoff:\t" + dist);	
		System.out.println("Min. Seq. Sep.:\t" + (minss==-1?"none":minss));
		System.out.println("Max. Seq. Sep.:\t" + (maxss==-1?"none":maxss));
		try {
			Model mod = new PdbFileModel(f, cc, ct, dist, minss, maxss);
			if(secondModel) {
				if(cmPane.addSecondModel(mod) == false) {
					showDifferentSizeContactMapsError();
				}else {
					compareStatus = true;
					cmPane.toggleCompareMode(compareStatus);
					Start.getPyMolAdaptor().loadStructure(mod.getTempPdbFileName(), mod.getPDBCode(), mod.getChainCode(), true);
					Start.getPyMolAdaptor().alignStructure(cmPane.getFirstModel().getPDBCode(), cmPane.getFirstModel().getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
				
				}
			} else {
				this.spawnNewViewWindow(mod);
			}
		} catch(ModelConstructionError e) {
			showLoadError(e.getMessage());
		}
	}	

	private void handleLoadFromCmFile(boolean secondModel) {
		
		if (secondModel == SECOND_MODEL && mod == null){
			this.showNoContactMapWarning();
		} else{
			LoadDialog dialog = new LoadDialog(this, "Load from Contact map file", new LoadAction(secondModel) {
				public void doit(Object o, String f, String ac, String cc, String ct, double dist, int minss, int maxss, String db, int gid) {
					View view = (View) o;
					view.doLoadFromCmFile(f, secondModel);
				}
			}, "", null, null, null, null, null, null, null, null);
			dialog.showIt();		  
		}
	}

	public void doLoadFromCmFile(String f, boolean secondModel) {
		System.out.println("Loading from contact map file "+f);
		try {
			Model mod = new ContactMapFileModel(f);
			if(secondModel) {
				if(cmPane.addSecondModel(mod) == false) {
					showDifferentSizeContactMapsError();
				}else {
					compareStatus = true;
					cmPane.toggleCompareMode(compareStatus);
					Start.getPyMolAdaptor().loadStructure(mod.getTempPdbFileName(), mod.getPDBCode(), mod.getChainCode(), true);
					Start.getPyMolAdaptor().alignStructure(cmPane.getFirstModel().getPDBCode(), cmPane.getFirstModel().getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
				
				}
			} else {
				this.spawnNewViewWindow(mod);
			}
		} catch(ModelConstructionError e) {
			showLoadError(e.getMessage());
		}
	}

	private void handleLoadFromFtp(boolean secondModel) {
		LoadDialog dialog = new LoadDialog(this, "Load from online PDB", new LoadAction(secondModel) {
			public void doit(Object o, String f, String ac, String cc, String ct, double dist, int minss, int maxss, String db, int gid) {
				View view = (View) o;
				view.doLoadFromFtp(ac, cc, ct, dist, minss, maxss, secondModel);
			}
		}, null, "", "", Start.DEFAULT_CONTACT_TYPE, String.valueOf(Start.DEFAULT_DISTANCE_CUTOFF), "", "", null, null);
		dialog.showIt();
	}

	public void doLoadFromFtp(String ac, String cc, String ct, double dist, int minss, int maxss, boolean secondModel) {
		System.out.println("Loading from online PDB");
		System.out.println("PDB code:\t" + ac);
		System.out.println("Chain code:\t" + cc);
		System.out.println("Contact type:\t" + ct);
		System.out.println("Dist. cutoff:\t" + dist);	
		System.out.println("Min. Seq. Sep.:\t" + (minss==-1?"none":minss));
		System.out.println("Max. Seq. Sep.:\t" + (maxss==-1?"none":maxss));	
		try{
			Model mod = new PdbFtpModel(ac, cc, ct, dist, minss, maxss);
			if(secondModel) {
				if(cmPane.addSecondModel(mod) == false) {
					showDifferentSizeContactMapsError();
				} else {
					compareStatus = true;
					cmPane.toggleCompareMode(compareStatus);
					Start.getPyMolAdaptor().loadStructure(mod.getTempPdbFileName(), mod.getPDBCode(), mod.getChainCode(), true);
					Start.getPyMolAdaptor().alignStructure(cmPane.getFirstModel().getPDBCode(), cmPane.getFirstModel().getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
					
				}
			 
			} else {
				this.spawnNewViewWindow(mod);
			}
		} catch(ModelConstructionError e) {
			showLoadError(e.getMessage());
		}
	}


	private void handleSaveToGraphDb() {
		if(this.mod == null) {
			showNoContactMapWarning();
		} else if(!Start.isDatabaseConnectionAvailable()) {
			showNoDatabaseConnectionWarning();
		} else {
			String ret = JOptionPane.showInputDialog(this, (Object) "<html>Name of database to save to:<br></html>");
			if(ret != null) {
				if(ret.length() > 0) {
					doSaveToGraphDb(ret);
				} else {
					JOptionPane.showMessageDialog(this, "No database name given. Contact map not saved.", "Warning", JOptionPane.WARNING_MESSAGE);
				}
			}
		}
	}

	public void doSaveToGraphDb(String dbName) {
		try {
			mod.writeToGraphDb(dbName);
			System.out.println("Saving contact map to database " + dbName + ".");
		} catch (SQLException e) {
			System.err.println("Error when trying to write to database " + dbName + ": " + e.getMessage());
		}
	}

	private void handleSaveToCmFile() {
		if(this.mod == null) {
			//System.out.println("No contact map loaded yet.");
			showNoContactMapWarning();
		} else {
			int ret = Start.getFileChooser().showSaveDialog(this);
			if(ret == JFileChooser.APPROVE_OPTION) {
				File chosenFile = Start.getFileChooser().getSelectedFile();
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
			int ret = Start.getFileChooser().showSaveDialog(this);
			if(ret == JFileChooser.APPROVE_OPTION) {
				File chosenFile = Start.getFileChooser().getSelectedFile();

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
			+ "Min Seq Sep: " + (mod.getMinSequenceSeparation()<1?"none":mod.getMinSequenceSeparation()) + "\n"
			+ "Max Seq Sep: " + (mod.getMaxSequenceSeparation()<1?"none":mod.getMaxSequenceSeparation()) + "\n"
			+ "\n"
			+ "Contact map size: " + mod.getMatrixSize() + "\n"
			+ "Unobserved Residues: " + mod.getNumberOfUnobservedResidues() + "\n"
			+ "Number of contacts: " + mod.getNumberOfContacts() + "\n"
			+ "Directed: " + (mod.isDirected()?"Yes":"No")
			+ "\n"
			+ "Sequence: " + s + "\n"
			+ "Secondary structure: " + mod.getSecondaryStructure().getComment();
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

	/* -------------------- View Menu -------------------- */
	
	/**
	 * 
	 */
	private void handleShowPdbResSers() {
		if(mod==null) {
			showNoContactMapWarning();
		} else if (!mod.has3DCoordinates()){
			showNo3DCoordsWarning();
		} else {
			showPdbSers = !showPdbSers;						// TODO: Move this to CMPane?
			cmPane.updateScreenBuffer();
			if(showPdbSers) {
				mmViewShowPdbResSers.setIcon(icon_selected);
			} else {
				mmViewShowPdbResSers.setIcon(icon_deselected);
			}
		}
	}
	
	/**
	 * 
	 */
	private void handleShowRulers() {
		if(mod==null) {
			showNoContactMapWarning();
		} else {
			showRulers = !showRulers;
			if(showRulers) {
				cmp.add(topRul, BorderLayout.NORTH);
				cmp.add(leftRul, BorderLayout.WEST);
				mmViewRulers.setIcon(icon_selected);
			} else {
				cmp.remove(topRul);
				cmp.remove(leftRul);
				mmViewRulers.setIcon(icon_deselected);
			}
			this.pack();
			this.repaint();
		}
	}
	
	/**
	 * 
	 */
	private void handleShowIconBar() {
		showIconBar = !showIconBar;
		if(showIconBar) {
			toolBar.setVisible(true);
			mmViewIconBar.setIcon(icon_selected);
		} else {
			toolBar.setVisible(false);
			mmViewIconBar.setIcon(icon_deselected);
		}
		this.repaint();		
	}	
	
	/**
	 * 
	 */
	private void handleShowNbhSizeMap() {
		if(mod==null) {
			showNoContactMapWarning();
		} else {
			showNbhSizeMap = !showNbhSizeMap;
			cmPane.toggleNbhSizeMap(showNbhSizeMap);
			if (showNbhSizeMap) {
				mmViewHighlightComNbh.setIcon(icon_selected);
			} else {
				mmViewHighlightComNbh.setIcon(icon_deselected);
			}
		}
	}

	/**
	 * 
	 */
	private void handleShowDensityMap() {
		if(mod==null) {
			showNoContactMapWarning();
		} else {
			showDensityMap = !showDensityMap;
			cmPane.toggleDensityMap(showDensityMap);
			if(showDensityMap) {
				mmViewShowDensity.setIcon(icon_selected);
			} else {
				mmViewShowDensity.setIcon(icon_deselected);
			}
		}
	}

	/**
	 * 
	 */
	private void handleShowDistanceMap() {
		if(mod==null) {
			showNoContactMapWarning();
		} else if (!mod.has3DCoordinates()){
			showNo3DCoordsWarning();
		} else {
			showDistanceMap = !showDistanceMap;
			cmPane.toggleDistanceMap(showDistanceMap);
			if(showDistanceMap) {
				mmViewShowDistMatrix.setIcon(icon_selected);
			} else {
				mmViewShowDistMatrix.setIcon(icon_deselected);
			}
		}
	}

	/* -------------------- Select menu -------------------- */
	
	private void handleSelectAll() {
		if(mod==null) {
			showNoContactMapWarning();
		} else {
			cmPane.selectAllContacts();
		}
	}
	
	private void handleSelectByResNum() {
		if(mod==null) {
			showNoContactMapWarning();
		} else {
			String selStr = (String)JOptionPane.showInputDialog(
					this,
					"Enter residue numbers to be selected:\n"
					+ "e.g. 1-5,7,10-11",
					"");
			if (selStr == null) {
				// user clicked cancel
				return;
			} else {
				if(!NodeSet.isValidSelectionString(selStr)) {
					showInvalidSelectionStringWarning();
					return;
				} else {
					int selectedContacts = cmPane.selectByResNum(selStr);
					System.out.println(selectedContacts + " contacts selected");					
				}
			}
		}
	}
	
	private void handleSelectHelixHelix() {
		if(mod==null) {
			showNoContactMapWarning();
		} else
		if(!mod.hasSecondaryStructure()) {
			showNoSecondaryStructureWarning();
		} else {
			cmPane.selectHelixHelix();
		}
	}
	
	private void handleSelectBetaBeta() {
		if(mod==null) {
			showNoContactMapWarning();
		} else
		if(!mod.hasSecondaryStructure()) {
			showNoSecondaryStructureWarning();
		} else {
			cmPane.selectBetaBeta();
		}
	}
	
	private void handleSelectInterSsContacts() {
		if(mod==null) {
			showNoContactMapWarning();
		} else
		if(!mod.hasSecondaryStructure()) {
			showNoSecondaryStructureWarning();
		} else {
			cmPane.selectInterSsContacts();
		}
	}
	
	private void handleSelectIntraSsContacts() {
		if(mod==null) {
			showNoContactMapWarning();
		} else
		if(!mod.hasSecondaryStructure()) {
			showNoSecondaryStructureWarning();
		} else {
			cmPane.selectIntraSsContacts();
		}
	}


	/* -------------------- Color menu -------------------- */
	
	private void handleColorSelect() {	
		JDialog chooseDialog = JColorChooser.createDialog(this, "Choose color", true, Start.getColorChooser(), 
				new ActionListener() {
					public void actionPerformed(ActionEvent e){
							Color c = Start.getColorChooser().getColor();
							currentPaintingColor = c;
							System.out.println("Color changed to " + c);
					}
				} , null); // no handler for cancel button
		chooseDialog.setVisible(true);
	}

	private void handleColorPaint() {
		if(mod==null) {
			showNoContactMapWarning();
		} else if(cmPane.getSelContacts().size() == 0) {
			showNoContactsSelectedWarning();
		} else {
			cmPane.paintCurrentSelection(currentPaintingColor);
			cmPane.resetSelections();
			cmPane.repaint();
		}
	}

	private void handleColorReset() {
		if(mod==null) {
			showNoContactMapWarning();
		} else {
			cmPane.resetUserContactColors();
		}
	}

	/* -------------------- Action menu -------------------- */
	
	/**
	 * 
	 */
	private void handleShowSelContacts3D() {
		if(mod==null) {
			showNoContactMapWarning();
		} else if(!mod.has3DCoordinates()) {
			showNo3DCoordsWarning();
		} else if(!Start.isPyMolConnectionAvailable()) {
			showNoPyMolConnectionWarning();
		} else if(cmPane.getSelContacts().size() == 0) {
			showNoContactsSelectedWarning();
		} else if (cmPane.hasSecondModel()){


			String firstModelContactColor = "color magenta, ";
			String secondModelContactColor = "color green, ";

			common = this.getSelCommonContactsInComparedMode();
			firstS = this.getSelFirstStrucInComparedMode();
			secondS = this.getSelSecondStrucInComparedMode();


			// only second structure contacts
			if (common == false && firstS == false && secondS == true){
				EdgeSet[] array = cmPane.getSelectedContacts();
				EdgeSet trueGreen = array[0];	// red contacts
				String selectionType = cmPane.getSecondModel().getChainCode();

				//disable all old objects and selections and enable the two actual objects 
				Start.getPyMolAdaptor().setView(mod.getPDBCode(), mod.getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
				
				// present contacts in second structure
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), "True"+selectionType, secondModelContactColor, pymolSelSerial, trueGreen, false, false);

				// unpresent contacts in main structure
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), selectionType, firstModelContactColor, pymolSelSerial, trueGreen, true, true);
				this.pymolSelSerial++;
			}

			// only first structure contacts
			else if (common == false && firstS == true && secondS == false){
				EdgeSet[] array = cmPane.getSelectedContacts();
				EdgeSet trueRed = array[0];	// red contacts
				String selectionType = mod.getChainCode();

				//disable all old objects and selections and enable the two actual objects 
				Start.getPyMolAdaptor().setView(mod.getPDBCode(), mod.getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
				
				// present contacts in main structure
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), "True"+selectionType, firstModelContactColor, pymolSelSerial, trueRed, false, false);
				// unpresent contacts in second structure
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), selectionType, secondModelContactColor, pymolSelSerial, trueRed, true, true);
				this.pymolSelSerial++;
			}

			// only first and second structure contacts
			else if (common == false && firstS == true && secondS == true){
				EdgeSet[] array = cmPane.getSelectedContacts();
				EdgeSet trueRed = array[0];		// red contacts
				EdgeSet trueGreen = array[1];	// green contacts
				String selectionType = mod.getChainCode() + cmPane.getSecondModel().getChainCode();

				// all contacts are either red or green. 
				// red contacts are n PyMol: red and not dashed in main structure && green and dashed in second structure
				// green contacts: analogous w.r.t. second structure and main structure

				//disable all old objects and selections and enable the two actual objects 
				Start.getPyMolAdaptor().setView(mod.getPDBCode(), mod.getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
				
				//present and unpresent contacts in main structure
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), "True"+selectionType, firstModelContactColor, pymolSelSerial, trueRed, false, false);	
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), selectionType, firstModelContactColor, pymolSelSerial,trueGreen, true, false);	

				// present and unpresent contacts in second structure
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), "True"+selectionType, secondModelContactColor, pymolSelSerial, trueGreen, false, false);
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), selectionType, secondModelContactColor, pymolSelSerial, trueRed, true, true);
				this.pymolSelSerial++;
			}

			// only common toggle mode			
			else if (common == true && firstS == false && secondS == false){
				EdgeSet[] array = cmPane.getSelectedContacts();
				EdgeSet trueRed = array[0];		// red contacts
				EdgeSet trueGreen = array[1];	// green contacts
				String selectionType = mod.getChainCode() + cmPane.getSecondModel().getChainCode();

				// no unpresent contacts

				//disable all old objects and selections and enable the two actual objects 
				Start.getPyMolAdaptor().setView(mod.getPDBCode(), mod.getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
				
				// present contacts in main and second structure
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), "True"+selectionType, firstModelContactColor, pymolSelSerial, trueRed, false, false);	
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), "True"+selectionType, secondModelContactColor, pymolSelSerial, trueGreen, false, true);	
				this.pymolSelSerial++;
			}

			// common and first structure mode == complete first structure
			else if (common == true && firstS == true && secondS == false){
				EdgeSet[] array = cmPane.getSelectedContacts();
				EdgeSet trueRed = array[0];		// red contacts
				EdgeSet trueGreen = array[1];	// green contacts
				String selectionType = mod.getChainCode();

				//disable all old objects and selections and enable the two actual objects 
				Start.getPyMolAdaptor().setView(mod.getPDBCode(), mod.getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
				
				//present contacts in main structure
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), "True"+selectionType, firstModelContactColor, pymolSelSerial, trueRed, false, false);	
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), "True"+selectionType, firstModelContactColor, pymolSelSerial, trueGreen, false, false);	

				// unpresent contacts in second structure
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), selectionType, secondModelContactColor, pymolSelSerial, trueGreen, true, false);
				//present contacts in second structure
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), "True"+selectionType, secondModelContactColor, pymolSelSerial, trueRed, false, true);
				this.pymolSelSerial++;
			}

			// common and second structure mode == complete second structure
			else if (common == true && firstS == false && secondS == true){
				EdgeSet[] array = cmPane.getSelectedContacts();
				EdgeSet trueRed = array[0];		// red contacts
				EdgeSet trueGreen = array[1];	// green contacts
				String selectionType =  cmPane.getSecondModel().getChainCode();

				//disable all old objects and selections and enable the two actual objects 
				Start.getPyMolAdaptor().setView(mod.getPDBCode(), mod.getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
				
				// unpresent contacts in main structure
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), selectionType, firstModelContactColor, pymolSelSerial, trueRed, true, false);	
				// present contacts in main structure
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), "True"+selectionType, firstModelContactColor, pymolSelSerial, trueGreen, false, false);	

				// present contacts in second structure
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), "True"+selectionType, secondModelContactColor, pymolSelSerial, trueGreen, false, false);
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), "True"+selectionType, secondModelContactColor, pymolSelSerial, trueRed, false, true);
				this.pymolSelSerial++;
			}
			else if (common == true && firstS == true && secondS == true){
				EdgeSet[] array = cmPane.getSelectedContacts();
				EdgeSet trueRed = array[0];		// red contacts
				EdgeSet trueGreen = array[1];	// green contacts
				EdgeSet trueCommon = array[2]; 	// common contacts
				String selectionType =  mod.getChainCode()+cmPane.getSecondModel().getChainCode();
				
				//disable all old objects and selections and enable the two actual objects 
				Start.getPyMolAdaptor().setView(mod.getPDBCode(), mod.getChainCode(), cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode());
				
				// present contacts in both structures
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), "True"+selectionType, firstModelContactColor, pymolSelSerial, trueCommon, false, false);	
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), "True"+selectionType, secondModelContactColor, pymolSelSerial, trueCommon, false, false);

				// present and unpresent contacts only in main structure
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), "True"+selectionType, firstModelContactColor, pymolSelSerial, trueRed, false, false);	
				Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), selectionType, firstModelContactColor, pymolSelSerial, trueGreen, true, false);	

				// present and unpresent contacts only in second struture
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), "True"+selectionType, secondModelContactColor, pymolSelSerial, trueGreen, false, false);
				Start.getPyMolAdaptor().edgeSelection(cmPane.getSecondModel().getPDBCode(), cmPane.getSecondModel().getChainCode(), selectionType, secondModelContactColor, pymolSelSerial, trueRed, true, true);
				this.pymolSelSerial++;
			}
		}
		else {
			EdgeSet[] array = cmPane.getSelectedContacts();
			EdgeSet contacts = array[0];
			String selectionType = mod.getChainCode();
			
			//disable all old objects and selections
			Start.getPyMolAdaptor().sendCommand("disable all");
			Start.getPyMolAdaptor().sendCommand("enable " + mod.getPDBCode()+mod.getChainCode());
			
			Start.getPyMolAdaptor().edgeSelection(mod.getPDBCode(), mod.getChainCode(), selectionType, "color magenta, ", pymolSelSerial, contacts, false, true);
			this.pymolSelSerial++;
		}
	}
		

	/**
	 * 
	 */
	private void handleShowDistance3D() {
		if(mod==null) {
			showNoContactMapWarning();
		} else if(!mod.has3DCoordinates()) {
			showNo3DCoordsWarning();
		} else if(!Start.isPyMolConnectionAvailable()) {				
			showNoPyMolConnectionWarning();
		} else {
			Start.getPyMolAdaptor().sendSingleEdge(mod.getPDBCode(), mod.getChainCode(), pymolSelSerial, cmPane.getRightClickCont());
			this.pymolSelSerial++;
		}
	}

	/**
	 * 
	 */
	private void handleShowTriangles3D() {
		if(mod==null) {
			showNoContactMapWarning();
		} else if(!mod.has3DCoordinates()) {
			showNo3DCoordsWarning();
		} else if(!Start.isPyMolConnectionAvailable()) {				
			showNoPyMolConnectionWarning();
		} else if(cmPane.getCommonNbh() == null) {
			showNoCommonNbhSelectedWarning();
		} else {
			Start.getPyMolAdaptor().showTriangles(mod.getPDBCode(), mod.getChainCode(), cmPane.getCommonNbh(), pymolNbhSerial); // TODO: get rid of NbhSerial
			this.pymolNbhSerial++;					
		}
	}

	/**
	 * 
	 */
	private void handleDeleteSelContacts() {
		if(mod==null) {
			showNoContactMapWarning();
		} else if(cmPane.getSelContacts().size() == 0) {
			showNoContactsSelectedWarning();
		} else {
			cmPane.deleteSelectedContacts();
		}
	}
	
	
	/* -------------------- Compare menu -------------------- */
	
	
	private void handleSelContactsInComparedMode(){		
		if(mod==null) {
			showNoContactMapWarning();
		} else
		if(!cmPane.hasSecondModel()) {
			showNoSecondContactMapWarning();
		} else {
			selCommonContactsInComparedMode = !selCommonContactsInComparedMode;
			cmPane.toggleCompareMode(selCommonContactsInComparedMode);
			if(selCommonContactsInComparedMode) {
				mmSelCommonContactsInComparedMode.setIcon(icon_selected);
			} else {
				mmSelCommonContactsInComparedMode.setIcon(icon_deselected);
			}
		}
		tbShowCommon.setSelected(selCommonContactsInComparedMode);
	}
	
	private void handleSelFirstStrucInComparedMode(){
		if(mod==null) {
			showNoContactMapWarning();
		} else
		if(!cmPane.hasSecondModel()) {
			showNoSecondContactMapWarning();
		} else {
			selFirstStrucInComparedMode = !selFirstStrucInComparedMode;
			cmPane.toggleCompareMode(selFirstStrucInComparedMode);
			tbShowFirstStructure.setSelected(selFirstStrucInComparedMode);
			if(selFirstStrucInComparedMode) {
				mmSelFirstStrucInComparedMode.setIcon(icon_selected);
			} else {
				mmSelFirstStrucInComparedMode.setIcon(icon_deselected);
			}
		}
		tbShowFirstStructure.setSelected(selFirstStrucInComparedMode);
	}
	
	private void handleSelSecondStrucInComparedMode(){
		if(mod==null) {
			showNoContactMapWarning();
		} else
		if(!cmPane.hasSecondModel()) {
			showNoSecondContactMapWarning();
		} else {
			selSecondStrucInComparedMode = !selSecondStrucInComparedMode;
			cmPane.toggleCompareMode(selSecondStrucInComparedMode);
			if(selSecondStrucInComparedMode) {
				mmSelSecondStrucInComparedMode.setIcon(icon_selected);
			} else {
				mmSelSecondStrucInComparedMode.setIcon(icon_deselected);
			}
		}
		tbShowSecondStructure.setSelected(selSecondStrucInComparedMode);
	}
	
	private void handleToggleDiffDistMap() {
		if(mod==null) {
			showNoContactMapWarning();
		} else if (!mod.has3DCoordinates()){
			showNo3DCoordsWarning();
		} else if (!cmPane.hasSecondModel()) {
			showNoSecondContactMapWarning();
		} else if (!cmPane.mod2.has3DCoordinates()) {
			showNo3DCoordsSecondModelWarning();
		} else {
			showDiffDistMap = !showDiffDistMap;
			cmPane.toggleDiffDistMap(showDiffDistMap);
			if(showDiffDistMap) {
				mmToggleDiffDistMap.setIcon(icon_selected);
			} else {
				mmToggleDiffDistMap.setIcon(icon_deselected);
			}
		}
	}
	
	/* -------------------- Help menu -------------------- */
	
	private void handleHelpWriteConfig() {
		try {
			Start.writeExampleConfigFile(Start.CONFIG_FILE_NAME);
			System.out.println("Writing example configuration file " + new File(Start.CONFIG_FILE_NAME).getAbsolutePath());
		} catch(IOException e) {
			System.err.println("Could not write configuration file " + new File(Start.CONFIG_FILE_NAME).getAbsolutePath());
		}
	}
	
	private void handleHelpAbout() {
		JOptionPane.showMessageDialog(this,
				"<html><center>" +
				"<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" +
				"Contact Map Viewer v"+Start.VERSION +
				"&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</p>" + 
				"<p>(C) AG Lappe 2007</p>" +
				"</center></html>",
				"About",
				JOptionPane.PLAIN_MESSAGE);
	}
	
	private void handleHelpHelp() {
		JOptionPane.showMessageDialog(this,
				"<html>" +
				"General<br>" +
				"- Click right mouse button in contact map for a context menu of available actions<br>" +
				"<br>" +
				"Square selection mode<br>" +
				"- Click on a contact to select it<br>" +
				"- Drag the mouse to select a rectangular area of contacts<br>" +
				"- Hold 'Ctrl' while selecting to add to the current selection<br>" +
				"- Click on a non-contact to reset the current selection<br>" +
				"<br>" +
				"Fill selection mode<br>" +
				"- Click on a contact to start a fill selection from that contact<br>" +
				"- Hold 'Ctrl' while selecting to add to the current selection<br>" +
				"<br>" +
				"Diagonal selection mode<br>" +
				"- Click to select all contacts along a diagonal<br>" +
				"- Click and drag to select multiple diagonals<br>" +
				"- Hold 'Ctrl' while selecting to add to the current selection<br>" +
				"<br>" +				
				"Node neighbourhood selection mode<br>" +
				"- Click on a residue in the ruler or in the diagonal to select its contacts<br>" +
				"- Click on a cell in the upper half to select all contacts of that pair of residues<br>" +
				"- Hold 'Ctrl' while selecting to add to the current selection<br>" +				
				"<br>" +
				"Show selected contacts in PyMol<br>" +
				"- Shows the currently selected contacts as edges in PyMol<br>" +
				"<br>" +
				"Show common neigbours<br>" +
				"- Click on a contact or non-contact to see the common neighbours for that pair of residues<br>" +
				"<br>" +
				"Show common neighbours in PyMol<br>" +
				"- Shows the last shown common neighbours as triangles in PyMol<br>" +
				"<br>" +
				"Delete selected contacts<br>" +
				"- Permanently deletes the selected contacts from the contact map<br>" +
				"</html>",
				"Help",
				JOptionPane.PLAIN_MESSAGE);
	}	
	
	/* -------------------- Warnings -------------------- */
	
	/** Shows a window with a warning message that no contact map is loaded yet */
	private void showNoContactMapWarning() {
		JOptionPane.showMessageDialog(this, "No contact map loaded yet", "Warning", JOptionPane.INFORMATION_MESSAGE);
	}
	
	/** Shows a window with a warning message that no secondary structure information is available */
	private void showNoSecondaryStructureWarning() {
		JOptionPane.showMessageDialog(this, "No secondary structure information available", "Warning", JOptionPane.INFORMATION_MESSAGE);
	}	
	
	/** Shows a window with a warning message that no second contact map is loaded yet (for comparison functions) */
	private void showNoSecondContactMapWarning() {
		JOptionPane.showMessageDialog(this, "No second contact map loaded yet", "Warning", JOptionPane.INFORMATION_MESSAGE);
	}	
	
	/** Warning dialog to be shown if a function is being called which requires 3D coordinates and they are missing */
	private void showNo3DCoordsWarning(){
		JOptionPane.showMessageDialog(this, "No 3D coordinates are associated with this contact map", "Warning", JOptionPane.INFORMATION_MESSAGE);
	}
	
	/** Warning dialog to be shown if a comparison function is being called which requires 3D coordinates and they are missing */
	private void showNo3DCoordsSecondModelWarning(){
		JOptionPane.showMessageDialog(this, "No 3D coordinates are associated with this contact map", "Warning", JOptionPane.INFORMATION_MESSAGE);
	}
	
	/** Error dialog to be shown if loading a model failed. */
	private void showLoadError(String message) {
		JOptionPane.showMessageDialog(this, "<html>Failed to load structure:<br>" + message + "</html>", "Load Error", JOptionPane.ERROR_MESSAGE);
	}
	
	/** Error dialog to be shown when trying to do a db operation without a db connection */
	private void showNoDatabaseConnectionWarning() {
		JOptionPane.showMessageDialog(this, "Failed to perform operation. No database connection available.", "Warning", JOptionPane.INFORMATION_MESSAGE);
	}
	
	/** Error dialog to be shown when trying to do a db operation without a db connection */
	private void showNoPyMolConnectionWarning() {
		JOptionPane.showMessageDialog(this, "Failed to perform operation. No Communication with PyMol available", "Warning", JOptionPane.INFORMATION_MESSAGE);
	}
	
	private void showNoCommonNbhSelectedWarning() {
		JOptionPane.showMessageDialog(this, "No common neighbourhood selected (Use Show Common Neighbours Mode)", "Warning", JOptionPane.INFORMATION_MESSAGE);		
	}
	
	private void showNoContactsSelectedWarning() {
		JOptionPane.showMessageDialog(this, "No contacts selected", "Warning", JOptionPane.INFORMATION_MESSAGE);				
	}
	
	private void showDifferentSizeContactMapsError() {
		JOptionPane.showMessageDialog(this, "The two contact maps can not be compared because they have different size.", "Can not compare contact maps", JOptionPane.ERROR_MESSAGE);		
	}
	
	private void showInvalidSelectionStringWarning() {
		JOptionPane.showMessageDialog(this, "Invalid selection string", "Warning", JOptionPane.INFORMATION_MESSAGE);				
	}

	/*---------------------------- public methods ---------------------------*/
	
	/** 
	 * Create and show a new view window, showing a contact map based on the given model and dispose the current one.
	 * TODO: Currently being abused as a constructor for new windows (something the real constructor should do).
	 */
	public void spawnNewViewWindow(Model mod) {

		String wintitle = "Contact Map of " + mod.getPDBCode() + " " + mod.getChainCode();
		View view = new View(mod, wintitle);
		if(view == null) {
			System.err.println("Error: Couldn't initialize contact map window");
			//System.exit(-1);
		}
		System.out.println("Contact map " + mod.getPDBCode() + " " + mod.getChainCode() + " loaded.");
		if(ContactMapPane.BG_PRELOADING) {
			view.cmPane.preloadBackgroundMaps();
		}

		if (Start.isPyMolConnectionAvailable() && mod.has3DCoordinates()) {
			// load structure in PyMol
			Start.getPyMolAdaptor().loadStructure(mod.getTempPdbFileName(), mod.getPDBCode(), mod.getChainCode(), false);
	
		}		
		// if previous window was empty (not showing a contact map) dispose it
		if(this.mod == null) {
			this.setVisible(false);
			this.dispose();
		}
	}
	
	/**
	 * Sets the comparison mode to the given state. If comparison mode is enabled, buttons
	 * for comparing structures become active otherwise they become inactive.
	 * @param state
	 */
	public void setComparisonMode(boolean state) {
		comparisonMode = !comparisonMode;
		if(comparisonMode) {
			tbShowCommon.setVisible(true);
			tbShowFirstStructure.setVisible(true);
			tbShowSecondStructure.setVisible(true);
		} else {
			tbShowCommon.setEnabled(false);
			tbShowFirstStructure.setVisible(false);
			tbShowSecondStructure.setVisible(false);			
		}
	}
	
	/* -------------------- getter methods -------------------- */
	// TODO: Make all these parameters of calls to CMPane
	/**
	 * Returns the current selection mode. The selection mode defines the meaning of user mouse actions.
	 * @return The current selection mode
	 */
	public int getCurrentSelectionMode(){
		return currentSelectionMode;
	}
	
	public boolean getCompareStatus(){
		return compareStatus;
	}
	
	/** 
	 * Returns whether showing the pdb residue serials is switched on 
	 */
	public boolean getShowPdbSers() {
		return showPdbSers;
	}

	/**  
	 * Returns whether showing the common neighbourhood size map is switched on
	 */
	public boolean getShowNbhSizeMap() {
		return showNbhSizeMap;
	}
	
	/** 
	 * Returns whether showing the density map is switched on
	 */
	public boolean getShowDensityMap() {
		return this.showDensityMap;
	}

	/** 
	 * Returns whether showing the distance map is switched on 
	 */
	public boolean getShowDistMap() {
		return this.showDistanceMap;
	}
	
	/** 
	 * Returns whether selection of contacts of both structures on compared contact map is switched on 
	 */
	public boolean getSelCommonContactsInComparedMode() {
		return this.selCommonContactsInComparedMode;
	}
	
	/** 
	 * Returns whether selection of contacts of the first structure on compared contact map is switched on
	 */
	public boolean getSelFirstStrucInComparedMode() {
		return this.selFirstStrucInComparedMode;
	}
	
	/** 
	 * Returns whether selection of contacts of the second structure on compared contact map is switched on
	 */
	public boolean getSelSecondStrucInComparedMode() {
		return this.selSecondStrucInComparedMode;
	}
	
	/** 
	 * Returns whether showing the difference distance map (in comparison mode) is switched on 
	 */
	public boolean getShowDiffDistMap() {
		return this.showDiffDistMap;
	}	
}

