package cmview;

import java.awt.Color;

/**
 * The current GUI state of a View instance. This currently includes the current selection mode, compare mode and whether various display options
 * are switched on or off (distance map, pdb residue numbers, etc).
 */
public class GUIState {

	/*------------------------------ constants ------------------------------*/
	
	protected static final SelMode INITIAL_SEL_MODE = SelMode.RECT;
	
	/*--------------------------- type definitions --------------------------*/
	
	// the selection mode
	protected enum SelMode {RECT, FILL, NBH, COMNBH, DIAG, COLOR, TOGGLE};
	
	/*--------------------------- member variables --------------------------*/
	
	// current gui state
	private View view;					// the parent view
	private SelMode selectionMode;		// current selection mode, modify using setSelectionMode
	private boolean showPdbSers;		// whether showing pdb serials is switched on
	private boolean showAlignmentCoords;// whether showing alignment coordinates in info corner is switched on
	private boolean showRulers;			// whether showing residue rulers is switched on
	private boolean showIconBar;		// whether showing the icon bar is switched on
	private boolean showNbhSizeMap;		// whether showing the common neighbourhood size map is switched on 
	private boolean showDensityMap;		// whether showing the density map is switched on
	private boolean showDistanceMap;	// whether showing the distance map is switched on
	private boolean showResidueScoringMap; // whether one of the generic residue scoring function maps is switched on
	private boolean showDiffDistMap; 	// whether showing the difference distance map is switched on
	private boolean showTFFctMap;       // whether showing the transferfunction based map is switched on
	private boolean showDeltaRankMap; 	// whether showing the delta rank map is switched on
	private boolean showRealTimeContacts;   // whether real time contacts mode is on or not
	
	// Maps for the bottom-left contact map background
	
	private boolean showBottomNbhSizeMap;
	private boolean showBottomDensityMap;
	private boolean showBottomDistanceMap;
	private boolean showBottomDeltaRankMap;
	private boolean showBottomDiffDistMap;
	private boolean showBottomTFFctMap;
	private boolean showBottomResidueScoringMap;
	
	private String bottomResidueScoringFunctionName;
	private String residueScoringFunctionName;
	
	// compare mode
	private boolean compareMode;		// whether we are in pairwise comparison mode (i.e. a second structure has been loaded)
	private Color paintingColor;		// current color for coloring contacts selected by the user
	private boolean showCommon;			// when true, common contacts displayed in compare mode (and selections are only for common)
	private boolean showFirst; 			// when true, contacts unique to first structure displayed in compare mode (and selections are only for first)
	private boolean showSecond; 		// when true, contacts unique to second structure displayed in compare mode (and selections are only for second)
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Initializes the GUI state with default values.
	 */
	GUIState(View view) {
		this.view = view;
		this.selectionMode = INITIAL_SEL_MODE;
		this.showPdbSers = Start.SHOW_PDB_RES_NUMS;
		this.showAlignmentCoords=Start.SHOW_ALIGNMENT_COORDS;
		this.showRulers=Start.SHOW_RULERS;
		this.showIconBar=Start.SHOW_ICON_BAR;
		this.showNbhSizeMap = false;
		this.showDensityMap=false;
		this.showDistanceMap=false;
		this.compareMode = false;	
		this.paintingColor = Color.blue;
		this.showCommon= true;
		this.showFirst= true;
		this.showSecond= true;
		this.showDiffDistMap = false;
		this.showDeltaRankMap = false;
		this.showResidueScoringMap = false;
		this.showTFFctMap = false;
		
		this.showBottomNbhSizeMap = false;
		this.showBottomDensityMap = false;
		this.showBottomDistanceMap = false;
		this.showBottomDiffDistMap = false;
		this.showBottomDeltaRankMap = false;
		this.showBottomResidueScoringMap = false;
		this.showBottomTFFctMap = false;
		
		this.showRealTimeContacts = Start.SHOW_CONTACTS_IN_REALTIME;
	}
	
	/*---------------------------- public methods ---------------------------*/

	/*---------------- getters ---------------*/

	/**
	 * @return the selectionMode
	 */
	protected SelMode getSelectionMode() {
		return selectionMode;
	}

	/**
	 * @return the showPdbSers
	 */
	protected boolean getShowPdbSers() {
		return showPdbSers;
	}
	
	/**
	 * @return the showAlignmentCoords
	 */
	protected boolean getShowAlignmentCoords() {
		return showAlignmentCoords;
	}	

	/**
	 * @return the showRulers
	 */
	protected boolean getShowRulers() {
		return showRulers;
	}

	/**
	 * @return the showIconBar
	 */
	protected boolean getShowIconBar() {
		return showIconBar;
	}

	/**
	 * @return the showNbhSizeMap
	 */
	protected boolean getShowNbhSizeMap() {
		return showNbhSizeMap;
	}
	
	protected boolean getShowBottomNbhSizeMap() {
		return showBottomNbhSizeMap;
	}
	
	/**
	 * @return the showDensityMap
	 */
	protected boolean getShowDensityMap() {
		return showDensityMap;
	}
	
	protected boolean getShowBottomDensityMap() {
		return showBottomDensityMap;
	}
	
	protected boolean getShowDeltaRankMap() {
		return showDeltaRankMap;
	}
	
	protected boolean getShowBottomDeltaRankMap() {
		return showBottomDeltaRankMap;
	}
	
	/**
	 * @return the showDistanceMap
	 */
	protected boolean getShowDistanceMap() {
		return showDistanceMap;
	}
	
	protected boolean getShowBottomDistanceMap() {
		return showBottomDistanceMap;
	}
	
	/**
	 * @return indicates whether a generic residueScoringFunctionMap is shown in the first background 
	 */
	
	protected boolean getShowResidueScoringMap() {
		return showResidueScoringMap;
	}
	
	/**
	 * @return indicates whether a generic residueScoringFunctionMap is shown in the second background 
	 */
	
	protected boolean getShowBottomResidueScoringMap() {
		return showBottomResidueScoringMap;
	}
	
	/**
	 * @return the compareMode
	 */
	protected boolean getCompareMode() {
		return compareMode;
	}

	/**
	 * @return the paintingColor
	 */
	protected Color getPaintingColor() {
		return paintingColor;
	}

	/**
	 * @return the showCommon
	 */
	protected boolean getShowCommon() {
		return showCommon;
	}

	/**
	 * @return the showFirst
	 */
	protected boolean getShowFirst() {
		return showFirst;
	}

	/**
	 * @return the showSecond
	 */
	protected boolean getShowSecond() {
		return showSecond;
	}

	/**
	 * @return the showDiffDistMap
	 */
	protected boolean getShowDiffDistMap() {
		return showDiffDistMap;
	}
	
	protected boolean getShowBottomDiffDistMap() {
		return showBottomDiffDistMap;
	}
	
	/**
	 * @return the showTFFctMap
	 */
	protected boolean getShowTFFctMap() {
		return showTFFctMap;
	}
	
	protected boolean getShowBottomTFFctMap() {
		return showBottomTFFctMap;
	}
	
	protected boolean getShowRealTimeContacts() {
		return showRealTimeContacts;
	}
	/*---------------- setters ---------------*/
	
	/**
	 * Sets the current selection mode. This sets the internal state variable and changes some gui components.
	 * Use getSelectionMode to retrieve the current state.
	 */
	protected void setSelectionMode(SelMode mode) {
		switch(mode) {
		case RECT: view.tbSquareSel.setSelected(true); break;
		case FILL: view.tbFillSel.setSelected(true); break;
		case NBH: view.tbNbhSel.setSelected(true); break;
		case COMNBH: view.tbShowComNbh.setSelected(true); break;
		case DIAG: view.tbDiagSel.setSelected(true); break;
		case COLOR : view.tbSelModeColor.setSelected(true); break;
		case TOGGLE: view.tbToggleContacts.setSelected(true); break;
		default: System.err.println("Error in setSelectionMode. Unknown selection mode " + mode); return;
		}
		this.selectionMode = mode;
	}
	
	protected void toggleRealTimeContacts() {
		this.showRealTimeContacts = !this.showRealTimeContacts;
	}
	
	/**
	 * @param showPdbSers the showPdbSers to set
	 */
	protected void setShowPdbSers(boolean showPdbSers) {
		this.showPdbSers = showPdbSers;
	}

	/**
	 * @param showRulers the showRulers to set
	 */
	protected void setShowRulers(boolean showRulers) {
		this.showRulers = showRulers;
	}

	/**
	 * @param showIconBar the showIconBar to set
	 */
	protected void setShowIconBar(boolean showIconBar) {
		this.showIconBar = showIconBar;
	}

	/**
	 * @param showNbhSizeMap the showNbhSizeMap to set
	 */
	protected void setShowNbhSizeMap(boolean showNbhSizeMap) {
		this.showNbhSizeMap = showNbhSizeMap;
	}
	
	protected void setShowBottomNbhSizeMap(boolean b) {
		this.showBottomNbhSizeMap = b;
	}
	
	/**
	 * @param showDensityMap the showDensityMap to set
	 */
	protected void setShowDensityMap(boolean showDensityMap) {
		this.showDensityMap = showDensityMap;
	}
	
	protected void setShowBottomDensityMap(boolean b) {
		this.showBottomDensityMap = b;
	}
	
	/**
	 * @param showDistanceMap the showDistanceMap to set
	 */
	protected void setShowDistanceMap(boolean showDistanceMap) {
		this.showDistanceMap = showDistanceMap;
	}

	protected void setShowBottomDistanceMap(boolean b) {
		this.showBottomDistanceMap= b;
	}
	
	protected void setShowResidueScoringMap(boolean b) {
		this.showResidueScoringMap = b;
	}
	
	protected void setShowBottomResidueScoringMap(boolean b) {
		this.showBottomResidueScoringMap = b;
	}
	
	/**
	 * @param compareMode the compareMode to set
	 */
	protected void setCompareMode(boolean compareMode) {
		this.compareMode = compareMode;
	}

	/**
	 * @param paintingColor the paintingColor to set
	 */
	protected void setPaintingColor(Color paintingColor) {
		this.paintingColor = paintingColor;
	}

	/**
	 * @param showCommon the showCommon to set
	 */
	protected void setShowCommon(boolean showCommon) {
		this.showCommon = showCommon;
	}

	/**
	 * @param showFirst the showFirst to set
	 */
	protected void setShowFirst(boolean showFirst) {
		this.showFirst = showFirst;
	}

	/**
	 * @param showSecond the showSecond to set
	 */
	protected void setShowSecond(boolean showSecond) {
		this.showSecond = showSecond;
	}

	/**
	 * @param showDiffDistMap the showDiffDistMap to set
	 */
	protected void setShowDiffDistMap(boolean showDiffDistMap) {
		this.showDiffDistMap = showDiffDistMap;
	}
	
	protected void setShowBottomDiffDistMap(boolean showDiffDistMap) {
		this.showBottomDiffDistMap = showDiffDistMap;
	}
	
	/**
	 * @param showTFFctMap the showTFFctMap to set
	 */
	protected void setShowTFFctMap(boolean showTFFctMap) {
		this.showTFFctMap = showTFFctMap;
	}
	protected void setShowBottomTFFctMap(boolean showTFFctMap) {
		this.showBottomTFFctMap = showTFFctMap;
	}

	
	public void setShowDeltaRankMap(boolean b) {
		this.showDeltaRankMap= b;
		
	}

	public void setShowBottomDeltaRankMap(boolean b) {
		this.showBottomDeltaRankMap = b;
	}

	public boolean getShowBackground() {
		return showDeltaRankMap || showDensityMap || showDiffDistMap || showDistanceMap || showResidueScoringMap || showTFFctMap;
	}
	public boolean getShowBottomBackground() {
		return showBottomDeltaRankMap || showBottomDensityMap || showBottomDiffDistMap || showBottomDistanceMap || showBottomResidueScoringMap || showBottomTFFctMap;
	}

	public void setResidueScoringFunctionName(boolean bottom, String name) {
		if (bottom) {
			bottomResidueScoringFunctionName = name;
		} else {
			residueScoringFunctionName = name;
		}
	}
	
	public String getResidueScoringFunctionName(boolean bottom) {
		if (bottom) {
			return bottomResidueScoringFunctionName;
		} else {
			return residueScoringFunctionName;
		}
	}
	
}
