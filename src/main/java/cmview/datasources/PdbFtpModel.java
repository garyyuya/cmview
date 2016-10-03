package cmview.datasources;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.io.mmcif.MMcifParser;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;

import owl.core.structure.*;
import owl.core.structure.graphs.RIGEnsemble;
import owl.core.structure.graphs.RIGGeometry;
import owl.core.util.FileFormatException;


import cmview.Start;

/** 
 * A contact map data model based on a structure loaded from a CIF file downloaded from pdb's ftp
 */
public class PdbFtpModel extends Model {
    
	
	private File cifFile;		// the cached cif file downloaded from the PDB
	private CiffileParser parser;
		
	/**
	 * Overloaded constructor to load the data.
	 * @throws IOException 
	 * @throws FileFormatException 
	 */
	public PdbFtpModel(String pdbCode, String edgeType, double distCutoff, int minSeqSep, int maxSeqSep) throws IOException, ModelConstructionError {
		this.edgeType = edgeType; 
		this.distCutoff = distCutoff;
		this.minSeqSep = minSeqSep;
		this.maxSeqSep = maxSeqSep;
		try {
			this.parser = new CiffileParser(pdbCode, Start.PDB_FTP_URL);
		} catch (FileFormatException e) {
			throw new ModelConstructionError(e.getMessage());
		}
		this.cifFile = parser.getCifFile();
	}
	
	public PdbFtpModel(File cifFile, String edgeType, double distCutoff, int minSeqSep, int maxSeqSep) throws IOException, ModelConstructionError {
		this.cifFile = cifFile;
		this.edgeType = edgeType; 
		this.distCutoff = distCutoff;
		this.minSeqSep = minSeqSep;
		this.maxSeqSep = maxSeqSep;
		try {
			this.parser = new CiffileParser(cifFile);
		} catch (FileFormatException e) {
			throw new ModelConstructionError(e.getMessage());
		}
	}
	
	public PdbFtpModel(Model mod) {
	    super(mod);
	}
	
	public PdbFtpModel copy() {
	    return new PdbFtpModel(this);
	}

	public File getCifFile() {
		return cifFile;
	}

	/**
	 * Loads the chain corresponding to the passed chain code identifier.
	 * @param pdbChainCode  pdb chain code of the chain to be loaded
	 * @param modelSerial  a model serial
	 * @throws ModelConstructionError
	 */
	@Override
	public void load(String pdbChainCode, int modelSerial) throws ModelConstructionError {
		load(pdbChainCode, modelSerial, false);
	}
	
	/**
	 * Loads the chain corresponding to the passed chain code identifier.
	 * If loadEnsembleGraph is true, the graph in this model will be the average graph of the ensemble of all models
	 * instead of the graph of the specified model only. The PdbChain object will still correspond to the given model number.
	 * @param pdbChainCode  pdb chain code of the chain to be loaded
	 * @param modelSerial  a model serial
	 * @param loadEnsembleGraph whether to set the graph in this model to the (weighted) ensemble graph of all models
	 * @throws ModelConstructionError
	 */
	public void load(String pdbChainCode, int modelSerial, boolean loadEnsembleGraph) throws ModelConstructionError {
		// load CIF file from online pdb
		try {
			//PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile,modelSerial);
			
			

	        MMcifParser parser = new SimpleMMcifParser();

	        SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();

	        // The Consumer builds up the BioJava - structure object.
	        // you could also hook in your own and build up you own data model.          
	        parser.addMMcifConsumer(consumer);

	        BufferedReader br = new BufferedReader(new FileReader(cifFile));
	       
	        parser.parse(br);
	         

	        // now get the protein structure.
	        Structure fullpdb = consumer.getStructure();
	        
	        br.close();
	        
	        
			
			if(pdbChainCode == null) {
				//this.pdb = fullpdb.getFirstChain();
				//pdbChainCode = this.pdb.getPdbChainCode(); 
				this.pdb = fullpdb.getChains(modelSerial-1).get(0);
				
				
			//} else if(!fullpdb.containsPdbChainCode(pdbChainCode)) {
			} else if(!fullpdb.hasPdbChain(pdbChainCode)) {
				throw new ModelConstructionError("Chain '" + pdbChainCode + "' not found");
			} else {
				//this.pdb = fullpdb.getChain(pdbChainCode);	
				this.pdb = fullpdb.getPolyChainByPDB(pdbChainCode, modelSerial-1);
			}
			//TODO handle secondary structure later
			//this.secondaryStructure = pdb.getSecondaryStructure();	// in case, dssp is n/a, use ss from pdb
			super.checkAndAssignSecondaryStructure();				// if dssp is a/, recalculate ss
			if(loadEnsembleGraph == false || this.parser.getModels().length == 1) {
			
				
				this.graph = Utils.getRIGraph(pdb, edgeType, distCutoff, modelSerial);
			} else {
				RIGEnsemble e = new RIGEnsemble(edgeType, distCutoff);
				try {
					e.loadFromMultiModelFile(this.cifFile);
					this.graph = e.getAverageGraph();
					this.graph.setPdbCode(this.pdb.getStructure().getPDBCode());
					this.graph.setChainCode(pdbChainCode);
					this.setIsGraphWeighted(true);
				} catch (IOException e1) {
					throw new ModelConstructionError("Error loading ensemble graph: " + e1.getMessage());
				}
			}
			
			// this.graph and this.residues are now available
			//TODO 4Corinna compute graph geometry and hand it over to ContactView
			if(Start.USE_CGAP) {
				//TODO handle Geom later
				//this.graphGeom = new RIGGeometry(this.graph, this.pdb);
				System.out.println("PdbFtpModel   GraphGeometry loaded");
				//this.graphGeom.printGeom();
			}

			// assign a loadedGraphId to this model
			String name = this.graph.getPdbCode()+this.graph.getChainCode();
			if (this.graph.getPdbCode().equals(PdbAsymUnit.NO_PDB_CODE)) {
				name = DEFAULT_LOADEDGRAPHID;
			} 
			this.loadedGraphID = Start.setLoadedGraphID(name, this);

			super.writeTempPdbFile();
			
			super.filterContacts(minSeqSep, maxSeqSep);
			super.printWarnings(pdbChainCode);
			
		} catch (PdbLoadException e) {
			System.err.println("Failed to load structure.");
			throw new ModelConstructionError(e.getMessage());
		} catch (IOException e) {
			System.err.println("Failed to load structure.");
			throw new ModelConstructionError(e.getMessage());
		} catch (FileFormatException e) {
			System.err.println("Failed to load structure.");
			throw new ModelConstructionError(e.getMessage());
		}
	}
	
	/**
	 * Gets chain codes for all chains being present in the source.
	 * 
	 * @throws GetterError
	 */
	public String[] getChains() throws PdbLoadException {
		return parser.getChains();
	}

	/**
	 * Gets model indices for all models being present in the source.
	 * 
	 * @return array of model identifiers, null if such thing
	 * @throws GetterError
	 */
	public Integer[] getModels() throws PdbLoadException {
		return parser.getModels();
	}
}
