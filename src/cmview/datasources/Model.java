package cmview.datasources;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.File;
import java.sql.SQLException;
import java.util.Collections;
import java.util.HashMap;

import actionTools.GetterError;
import actionTools.TinkerStatusNotifier;

import proteinstructure.AAinfo;
import proteinstructure.Alignment;
import proteinstructure.IntPairSet;
import proteinstructure.Pdb;
import proteinstructure.PdbLoadError;
import proteinstructure.RIGCommonNbhood;
import proteinstructure.RIGEdge;
import proteinstructure.RIGNbhood;
import proteinstructure.RIGNode;
import proteinstructure.RIGraph;
import proteinstructure.SecondaryStructure;
import tinker.TinkerError;
import tinker.TinkerRunner;

import cmview.Start;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * A contact map data model. Derived classes have to implement the constructor
 * in which the structure is loaded, the member variables are set appropriately
 * and a temporary pdb file with the atom lines is written. All other methods
 * should be implemented here.
 */
public abstract class Model {

	public static final String DEFAULT_LOADEDGRAPHID = "CM";

	/* These members' only purpose is to temporarily store some values passed to the constructor and use them in load().
	   Eventually these should be parameters of load directly, but for the moment we keep them here in order to keep
	   the interface unchanged.
	   Now, they're also useful to reload a model with the same parameters after reconstruction
	*/
	public String edgeType;
	public double distCutoff;
	public int minSeqSep;
	public int maxSeqSep;
	
	
	/*--------------------------- member variables --------------------------*/

	// structure and contact map data
	protected Pdb pdb; // this can be null if there are no 3D coordinates
						// available
	protected RIGraph graph; // currently every model must have a valid graph
								// object
	protected HashMap<Pair<Integer>, Double> distMatrix; // the scaled [0-1]
															// distance matrix
															// TODO: Could move
															// this to
															// ContactMapPane

	private File tempPdbFile; // the file with the atomic coordinates to be
								// loaded into pymol

	protected String loadedGraphID; // the unique identifier of a user-loaded
									// graph, we assign this whenever the user

	// loads a graph/structure into the viewer, this is being used to identify
	// the
	// graph in an alignment object and in pymol

	/*----------------------------- constructors ----------------------------*/

	// see derived classes for main constructors
	// In implemented constructors, the variables pdb and graph have to be
	// initialized
	// and the following private methods have to be called in turn.

	// needed as implicit super constructor for derived classes
	public Model() {
	}

	/**
	 * Create a new model as a shallow copy of the given model
	 * 
	 * @param mod
	 *            the original model
	 */
	public Model(Model mod) {
		this.pdb = mod.pdb;
		this.graph = mod.graph;
		this.distMatrix = mod.distMatrix;
		this.tempPdbFile = mod.tempPdbFile;
		this.loadedGraphID = mod.loadedGraphID;
	}

	/**
	 * Creates a new model as a shallow copy of this model.
	 * 
	 * @return the new model
	 */
	public abstract Model copy();

	public abstract void load(String pdbChainCode, int modelSerial)
			throws ModelConstructionError;

	/*---------------------------- private methods --------------------------*/

	/**
	 * Write temporary PDB file with atom lines for the current structure.
	 * has3DCoordinates() must be true before calling this (i.e. pdb not null)
	 */
	protected void writeTempPdbFile() {
		try {
			pdb.dump2pdbfile(getTempPdbFileName());
		} catch (IOException e) {
			System.err.println("Error writing temporary PDB file "
					+ getTempPdbFileName());
		}
	}

	/**
	 * Filter out unwanted contacts and initializes the seqSep variable.
	 * Warning: Not tested for directed graphs
	 */
	protected void filterContacts(int minSeqSep, int maxSeqSep) {
		if (minSeqSep > 0) {
			this.graph.restrictContactsToMinRange(minSeqSep);
		}
		if (maxSeqSep > 0) {
			this.graph.restrictContactsToMaxRange(maxSeqSep);
		}
	}

	/** Print some warnings if necessary */
	protected void printWarnings(String oldChainCode) {
		if (getNumberOfUnobservedResidues() > 0)
			System.out.println("Warning: " + getNumberOfUnobservedResidues()
					+ " unobserved or non-standard residues");
		if (!oldChainCode.equals(getChainCode()))
			System.out.println("Warning: Chain " + oldChainCode
					+ " is now called " + getChainCode());
	}

	/**
	 * Check whether secondary structure needs to be assigned and run dssp if so
	 */
	protected void checkAndAssignSecondaryStructure() {
		if (has3DCoordinates()) { // otherwise we can't (re)assign secondary
									// structure
			if (Start.isDsspAvailable()) {
				System.out
						.println("(Re)assigning secondary structure using DSSP");
				try {
					pdb.runDssp(Start.DSSP_EXECUTABLE, Start.DSSP_PARAMETERS);
				} catch (IOException e) {
					System.err.println("Failed to assign secondary structure: "
							+ e.getMessage());
				}
			}
		}
	}

	/*---------------------------- public methods ---------------------------*/

	/**
	 * Gets chain codes for all chains being present in the source.
	 * 
	 * @throws GetterError
	 */
	public String[] getChains() throws PdbLoadError {
		return pdb.getChains();
	}

	/**
	 * Gets model indices for all models being present in the source.
	 * 
	 * @return array of model identifiers, null if such thing
	 * @throws GetterError
	 */
	public Integer[] getModels() throws PdbLoadError {
		return pdb.getModels();
	}

	/**
	 * Get pdb
	 * 
	 * @return pdb
	 */

	public Pdb getPdb() {
		return pdb;
	}

	/** Returns the size of the data matrix */
	public int getMatrixSize() {
		return this.graph.getFullLength();
	}

	/** Returns the contacts as an EdgeSet */
	public IntPairSet getContacts() {
		IntPairSet contacts = new IntPairSet();
		for (RIGEdge edge : this.graph.getEdges()) {
			Pair<RIGNode> pair = this.graph.getEndpoints(edge);
			contacts.add(new Pair<Integer>(pair.getFirst().getResidueSerial(),
					pair.getSecond().getResidueSerial()));
		}
		return contacts;
	}

	/** Returns the graph object */
	public RIGraph getGraph() {
		return this.graph;
	}

	/** Returns the number of contacts */
	public int getNumberOfContacts() {
		return graph.getEdgeCount();
	}

	/** Returns true if the graph is directed, false otherwise */
	public boolean isDirected() {
		return graph.isDirected();
	}

	/**
	 * Returns the unique identifier for this model
	 * 
	 * @return
	 */
	public String getLoadedGraphID() {
		return loadedGraphID;
	}

	/** Returns the pdb code of the underlying structure */
	public String getPDBCode() {
		return graph.getPdbCode();
	}

	/**
	 * Returns the internal chain code of the underlying structure. Note that
	 * the internal chain code may be different from the pdb chain code given
	 * when loading the structure. For a database model, the internal chain code
	 * is the one given by the database (Pdbase). When loading from pdb file,
	 * the internal chain code is the same as the pdb chain code or 'A' if the
	 * pdb chain code is empty.
	 */
	public String getChainCode() {
		return graph.getChainCode(); // gets the internal chain code (may be !=
										// pdb chain code)
	}

	/**
	 * Returns the pdb model number (e.g. for NMR structures) for the underlying
	 * graph.
	 * 
	 * @return the model number
	 */
	public int getPdbModelNumber() {
		return graph.getModel();
	}

	/**
	 * Returns the CASP target number, or 0 if graph is not a CASP graph
	 * 
	 * @return
	 */
	public int getTargetNum() {
		return graph.getTargetNum();
	}

	/**
	 * Returns the CASP group number, or 0 if graph is not a CASP graph
	 * 
	 * @return
	 */
	public int getGroupNum() {
		return graph.getGroupNum();
	}

	/**
	 * Returns the CASP model number, or 0 if graph is not a CASP graph
	 * 
	 * @return
	 */
	public int getCaspModelNum() {
		return graph.getCaspModelNum();
	}

	/**
	 * Returns the contact type
	 */
	public String getContactType() {
		return graph.getContactType();
	}

	/**
	 * Returns the distance cutoff
	 */
	public double getDistanceCutoff() {
		return graph.getCutoff();
	}

	/** Returns the sequence separation of the current graph */
	public int getMinSequenceSeparation() {
		return this.graph.getMinSeqSep();
	}

	/** Returns the sequence separation of the current graph */
	public int getMaxSequenceSeparation() {
		return this.graph.getMaxSeqSep();
	}

	/**
	 * Returns the sequence for the current structure or an empty string if no
	 * sequence information is available.
	 */
	public String getSequence() {
		return graph.getSequence();
	}

	/**
	 * Returns the temporary pdb-style file with atomic coordinates. When called
	 * for the first time, this method creates the File object and marks it to
	 * be deleted on exit.
	 */
	protected File getTempPdbFile() {
		if (tempPdbFile == null) {
			tempPdbFile = new File(Start.TEMP_DIR, getLoadedGraphID() + ".pdb");
			tempPdbFile.deleteOnExit(); // will delete the file when the VM is
										// closed. Unless it crashes
		}
		return tempPdbFile;
	}

	/** Returns the name of the temporary pdb file */
	public String getTempPdbFileName() {
		return getTempPdbFile().getAbsolutePath();
	}

	/**
	 * Returns the number of unobserved or non-standard residues in the
	 * structure. The contacts of these residues are ignored.
	 */
	public int getNumberOfUnobservedResidues() {
		return (graph.getFullLength() - graph.getObsLength());
	}

	/** Write the current contact map to a contact map file */
	public void writeToContactMapFile(String fileName) throws IOException {
		try {
			this.graph.write_graph_to_file(fileName);
		} catch (IOException e) {
			System.err.println("Error when trying to write contact map file");
			throw e;
		}
	}

	public void writeToCaspRRFile(String fileName) throws IOException {
		try {
			this.graph.writeToCaspRRFile(fileName);
		} catch (IOException e) {
			System.err.println("Error when trying to write CASP RR file");
			throw e;
		}
	}

	/** Write the current contact map to a graph database */
	public void writeToGraphDb(String dbName) throws SQLException {
		graph.write_graph_to_db(Start.getDbConnection(), dbName);
	}

	/**
	 * Returns the three letter residue type for the given residue serial. Will
	 * always return the residue type from the sequence of the RIGraph,
	 * independently of whether a RIGNode exists in the RIGraph for this residue
	 * serial.
	 * 
	 * @param resser
	 *            the residue serial
	 * @return A string with the three letter residue type of the residue with
	 *         serial resser
	 */
	public String getResType(int resser) {
		return AAinfo.oneletter2threeletter(String.valueOf(this.getSequence()
				.charAt(resser - 1)));

		// NOTE, we used to take the residue type from the RIGNode, but this
		// causes problems
		// when there are missing RIGNodes, which happens in 2 cases:
		// a) unobserved residues
		// b) in some contact types with no atoms for some residues (e.g. ALA in
		// Cg contact type:
		// there can't be any edges for ALA in Cg and thus also no RIGNode is
		// added to the RIGraph)
		// RIGNode node = graph.getNodeFromSerial(resser);
		// if (node==null) return null;
		// return node.getResidueType();
	}

	/**
	 * Returns the one letter residue type for the given residue serial.
	 * 
	 * @param resser
	 *            the residue serial
	 * @return A string with the one letter residue type of the residue with
	 *         serial resser
	 * @see {@#getResType(int)}
	 */
	public String getResType1Letter(int resser) {
		return String.valueOf(this.getSequence().charAt(resser - 1));
	}

	/**
	 * Returns the RIGNode for the given residue serial
	 * 
	 * @param resser
	 * @return
	 */
	public RIGNode getNodeFromSerial(int resser) {
		return this.graph.getNodeFromSerial(resser);
	}

	public RIGNbhood getNbhood(int i_resser) {
		return graph.getNbhood(graph.getNodeFromSerial(i_resser));
	}

	public RIGCommonNbhood getCommonNbhood(int i_resser, int j_resser) {
		return graph.getCommonNbhood(graph.getNodeFromSerial(i_resser), graph
				.getNodeFromSerial(j_resser));
	}

	/**
	 * To get pdb residue serial from internal residue serial May only be called
	 * if has3DCoordinates is true
	 * 
	 * @param resser
	 * @return
	 */
	public String getPdbResSerial(int resser) {
		return pdb.get_pdbresser_from_resser(resser);
	}

	public HashMap<Pair<Integer>, Integer> getAllCommonNbhSizes() {
		return graph.getAllCommonNbhSizes();
	}

	public boolean containsEdge(Pair<Integer> cont) {
		return graph.containsEdge(graph.findEdge(graph.getNodeFromSerial(cont
				.getFirst()), graph.getNodeFromSerial(cont.getSecond())));
	}

	public void addEdge(Pair<Integer> cont) {
		graph.addEdgeIJ(cont.getFirst(), cont.getSecond());

	}

	public void removeEdge(Pair<Integer> cont) {
		graph.removeEdge(this.graph.findEdge(graph.getNodeFromSerial(cont
				.getFirst()), graph.getNodeFromSerial(cont.getSecond())));
	}

	/**
	 * Initialises the distMatrix member. The matrix is scaled to doubles from
	 * zero to one. Returns the scaled value of the current distance cutoff. May
	 * only be called if has3DCoordinates is true.
	 */
	public double initDistMatrix() {
		HashMap<Pair<Integer>, Double> distMatrixRes = this.pdb
				.calculate_dist_matrix(Start.DIST_MAP_CONTACT_TYPE);
		double max = Collections.max(distMatrixRes.values());
		double min = Collections.min(distMatrixRes.values());
		distMatrix = new HashMap<Pair<Integer>, Double>(); // TODO: Use old
															// matrix to save
															// memory
		for (Pair<Integer> cont : distMatrixRes.keySet()) {
			distMatrix.put(cont, (distMatrixRes.get(cont) - min) / (max - min));
		}
		double dist = (graph.getCutoff() - min) / (max - min);
		return dist;
	}

	/**
	 * Returns the current distance matrix. Before initDistMatrix has been
	 * called this will be null.
	 * 
	 * @return A map assigning to each edge the corresponding distance (scaled
	 *         to [0;1]).
	 */
	public HashMap<Pair<Integer>, Double> getDistMatrix() {
		return distMatrix;
	}

	/**
	 * Returns the difference of the distance maps of this model and the given
	 * second model. All distances are based on C-alpha atoms.
	 * 
	 * @param ali
	 *            an alignment between this model's sequence and secondModel's
	 *            sequence
	 * @param secondModel
	 *            the second model to compare agains
	 * @return A map assigning to each edge the corresponding value in the
	 *         difference distance matrix or null on error.
	 */
	public HashMap<Pair<Integer>, Double> getDiffDistMatrix(Alignment ali,
			Model secondModel) {
		/*
		 * TODO: Also force c-alpha for simple distance maps? Throw proper
		 * exceptions instead of returning null? Use real matrix?
		 */
		double diff, min, max;
		if (!this.has3DCoordinates() || !secondModel.has3DCoordinates()) {
			System.err
					.println("Failed to compute difference distance map. No 3D coordinates.");
			return null; // distance doesn't make sense without 3D data
		}

		String name1 = this.getLoadedGraphID();
		String name2 = secondModel.getLoadedGraphID();
		HashMap<Pair<Integer>, Double> diffDistMatrix = this.pdb
				.getDiffDistMap(Start.DIST_MAP_CONTACT_TYPE, secondModel.pdb,
						Start.DIST_MAP_CONTACT_TYPE, ali, name1, name2);

		if (diffDistMatrix == null) {
			System.err.println("Failed to compute difference distance map.");
		} else {
			max = Collections.max(diffDistMatrix.values());
			min = Collections.min(diffDistMatrix.values());
			if (max == min)
				System.err
						.println("Failed to scale difference distance matrix. Matrix is empty or all matrix entries are the same.");
			else {
				// scale matrix to [0;1]
				for (Pair<Integer> e : diffDistMatrix.keySet()) {
					diff = diffDistMatrix.get(e);
					diffDistMatrix.put(e, (diff - min) / (max - min));
				}
			}
		}
		return diffDistMatrix;
	}

	/**
	 * Returns true if this model contains 3D coordinates. Certain other methods
	 * can only be used if 3D coordinates are available.
	 * 
	 * @return true if this model has 3D coordinates, false otherwise
	 */
	public boolean has3DCoordinates() {
		return (pdb != null);
	}

	/**
	 * Gets the pdb coordinates.
	 * 
	 * @return the coordinates, returns null if there is no such information
	 *         available.
	 * @see #has3DCoordinates()
	 */
	public Pdb get3DCoordinates() {
		return pdb;
	}

	/**
	 * Returns true if full sequence information is available
	 * 
	 * @return true if sequence information is available, false otherwise
	 */
	public boolean hasSequence() {
		return graph.hasSequence();
	}

	// secondary structure related methods

	/**
	 * Returns true if this model has secondary structure information.
	 * 
	 * @return true if secondary structure information is available, false
	 *         otherwise
	 */
	public boolean hasSecondaryStructure() {
		if (this.pdb == null) {
			return false;
		}
		return this.pdb.hasSecondaryStructure();
	}

	/**
	 * Returns the secondary structure annotation object of this model.
	 * 
	 * @return the secondary structure annotation object
	 */
	public SecondaryStructure getSecondaryStructure() {

		// We don't need to try to get the secondary structure from graph as at
		// the moment there
		// is no way we could have a graph without 3D coordinates that has
		// secondary structure
		// annotation (neither in our database nor in our file formats we have a
		// way to store the secondary structure annotation)

		// assert pdb.getSecondaryStructure()==graph.getSecondaryStructure();

		if (pdb != null && pdb.hasSecondaryStructure()) {
			return pdb.getSecondaryStructure();
		} else {
			return new SecondaryStructure("");
		}
	}

	// end of secondary structure related methods

	/**
	 * Calculates and returns the density matrix for the current model.
	 */
	public double[][] getDensityMatrix() {
		int size = getMatrixSize();
		double[][] d = new double[size][size]; // density matrix

		// initialize matrix with contacts
		for (RIGEdge cont : graph.getEdges()) {
			Pair<RIGNode> pair = graph.getEndpoints(cont);
			d[pair.getFirst().getResidueSerial() - 1][pair.getSecond()
					.getResidueSerial() - 1] = 1;
		}

		// fill diagonal - avoids artefacts near main diagonal
		for (int i = 0; i < size; i++) {
			d[i][i] = 1;
		}

		// starting from second diagonal, fill uppper matrix with contact counts
		// and calculate diagonal averages
		double[] avg = new double[size];
		double[] std = new double[size];
		for (int j = 2; j < size; j++) {
			double sum = 0;
			for (int i = 0; i < size - j; i++) {
				d[i][i + j] = d[i + 1][i + j] + d[i][i + j - 1]
						- d[i + 1][i + j - 1] + d[i][i + j];
				// down left left+down contact
				sum += d[i][i + j]; // sum of diagonal values
			}
			avg[j] = sum / (size - j); // diagonal average

			// standard deviation
			sum = 0;
			for (int i = 0; i < size - j; i++) {
				sum += Math.pow(d[i][i + j] - avg[j], 2);
			}
			std[j] = Math.sqrt(sum / (size - j));
			// System.out.println("avg(" + j + ") = " + avg[j]);
			// System.out.println("std(" + j + ") = " + std[j]);

		}

		// diagonal average and standard deviation for first off-diagonal
		double sum = 0;
		for (int i = 0; i < size - 1; i++) {
			sum = sum + d[i][i + 1];
		}
		avg[1] = sum / (size - 1);
		sum = 0;
		for (int i = 0; i < size - 1; i++) {
			sum += Math.pow(d[i][i + 1] - avg[1], 2);
		}
		std[1] = Math.sqrt(sum / (size - 1));

		// calculate z-score density
		for (int i = 0; i < size - 1; i++) {
			for (int j = i + 1; j < size; j++) {
				if (std[j - i] == 0) {
					d[i][j] = 0;
				} else {
					d[i][j] = (d[i][j] - avg[j - i]) / std[j - i];
				}
			}
		}

		// reset diagonal
		for (int i = 0; i < size; i++) {
			d[i][i] = 0;
		}

		// map values to [0,1]
		double max = 0;
		double min = d[0][0];
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (d[i][j] > max)
					max = d[i][j];
				if (d[i][j] < min)
					min = d[i][j];
			}
		}
		max = Math.min(max, 3); // cut off at 3
		min = Math.max(min, -3); // cut off at -3
		assert (max >= 0);
		assert (min <= 0);
		if (max - min > 0) {
			for (int i = 0; i < size; i++) {
				for (int j = i; j < size; j++) {
					d[i][j] = (d[i][j] - min) / (max - min);
				}
			}
		}
		return d;
	}

	public Pdb runTinker(TinkerStatusNotifier tinkerStatusNotifier, TinkerRunner.PARALLEL parallel,
			TinkerRunner.REFINEMENT refinement, int models, String tmpDir) {
		TinkerRunner tinker = null;
		try {
			tinker = new TinkerRunner(Start.TINKER_BINPATH,
					"/project/StruPPi/Software/tinker/amber99.prm");
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		try {
			tinker.setNotifier(tinkerStatusNotifier);
			if (Start.TINKER_TEMP_DIR != null) {
				tinker.setTmpDir(Start.TINKER_TEMP_DIR);	
			} else {
				tinker.setTmpDir(Start.TEMP_DIR);
			}
			return tinker.reconstruct(getSequence(), new RIGraph[] { graph },
					null, false, models, true);
		} catch (TinkerError e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}

	
}