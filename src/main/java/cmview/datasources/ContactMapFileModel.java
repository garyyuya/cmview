package cmview.datasources;
import cmview.Start;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmcif.MMcifParser;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;

import owl.core.structure.*;
import owl.core.structure.graphs.FileRIGraph;
import owl.core.util.FileFormatException;

/** 
 * A contact map data model based on a graph loaded from a Contact map file.
 */
public class ContactMapFileModel extends Model {
	
	/**
	 * Overloaded constructor to load the data from file.
	 * @throws ModelConstructionError 
	 */
	public ContactMapFileModel(String fileName) throws ModelConstructionError {
			
		// load Contact graph from file
		try {
			
			this.graph = new FileRIGraph(fileName);
			this.isGraphWeighted = graph.hasWeightedEdges();
			this.secondaryStructure = graph.getSecondaryStructure(); // take ss from graph
			
			String pdbCode = graph.getPdbCode();
			String pdbChainCode = graph.getPdbChainCode();
			int modelSerial = graph.getModel();
			
			// check whether sequence info exists
			if(graph.getSequence().equals("")) { //TODO this shouldn't happen since we don't allow blank sequences in CM files, but we keep as another check doesn't harm
				throw new ModelConstructionError("File contains no sequence information.");
			}
			
			// assign a loadedGraphId to this model
			String name = this.graph.getPdbCode()+this.graph.getChainCode();
			if (this.graph.getPdbCode().equals(PdbAsymUnit.NO_PDB_CODE)) {
				name = DEFAULT_LOADEDGRAPHID;
			} 
			this.loadedGraphID = Start.setLoadedGraphID(name, this);
			
			// load structure from pdbase/online if possible
			if(!pdbCode.equals(PdbAsymUnit.NO_PDB_CODE) && !pdbChainCode.equals(PdbChain.NO_PDB_CHAIN_CODE)) {
				// we try to load from online cif file
				try {
					File cifFile = new File(Start.TEMP_DIR,pdbCode + ".cif");
					PdbAsymUnit.grabCifFile(null, Start.PDB_FTP_URL, pdbCode, cifFile, true);
					//PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile, modelSerial);
					MMcifParser parser = new SimpleMMcifParser();

			        SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();

			        // The Consumer builds up the BioJava - structure object.
			        // you could also hook in your own and build up you own data model.          
			        parser.addMMcifConsumer(consumer);

			       
			        parser.parse(new BufferedReader(new FileReader(cifFile)));
			         

			        // now get the protein structure.
			        Structure fullpdb = consumer.getStructure();
					//this.pdb = fullpdb.getChain(pdbChainCode);
					this.pdb = fullpdb.getPolyChainByPDB(pdbChainCode, modelSerial-1);
					super.writeTempPdbFile(); // this doesn't make sense without a pdb object
					
					
					if(this.pdb != null) {
					this.secondaryStructure = Utils.convertSecondStruc(fullpdb, pdb);
					//this.secondaryStructure = pdb.getSecondaryStructure(); 
				}	
					
					
					
				}  	catch(IOException e) {
					System.err.println("Failed to load structure because of error while downloading/reading the CIF file: "+e.getMessage());
					pdb = null;
				} 
				// if pdb creation failed then pdb=null					

				
				
				// if structure is available, and has secondary structure annotation, use it	
			} else {
				System.out.println("No pdb code and/or chain code found. Can not load structure.");
			} 
			
			//super.filterContacts(seqSep);	// currently not allowed to filter contacts
			//super.printWarnings(chainCode); // doesn't make sense here
			
		} catch (IOException e) {
			System.err.println("Error while trying to load graph from contact map file.");
			throw new ModelConstructionError(e.getMessage());
		} catch (FileFormatException e){
			System.err.println("Error while trying to load graph from contact map file. Wrong graph file format.");
			throw new ModelConstructionError(e.getMessage());			
		}
		

	}
	
	public ContactMapFileModel(Model mod) {
	    super(mod);
	}
	
	public ContactMapFileModel copy() {
	    return new ContactMapFileModel(this);
	}

	/**
	 * The loading of the contact map is implemented in the constructor not in 
	 * this function. This function essentially does'nt do anything!
	 * @param pdbChainCode pdb chain code of the chain to be loaded (ignored!)
	 * @param modelSerial  a model serial
	 * @throws ModelConstructionError
	 */
	@Override
	public void load(String pdbChainCode, int modelSerial) throws ModelConstructionError {
		return;
	}
	
}
