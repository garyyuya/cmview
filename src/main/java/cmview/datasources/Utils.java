package cmview.datasources;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.TreeMap;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.contact.GroupContact;
import org.biojava.nbio.structure.contact.GroupContactSet;
import org.biojava.nbio.structure.secstruc.SecStrucCalc;
import org.biojava.nbio.structure.secstruc.SecStrucTools;

import edu.uci.ics.jung.graph.util.Pair;
import owl.core.sequence.Sequence;
import owl.core.structure.AaResidue;
import owl.core.structure.Atom;
import owl.core.structure.PdbChain;
import owl.core.structure.Residue;
import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.Interval;
import owl.core.util.IntervalSet;
import owl.core.util.MySQLConnection;


public class Utils {
	

	public static RIGraph getRIGraph(Chain chain, String edgeType, double cutoff, int modelSerial) {
		
		String[] atomNames = null;
		
		if (edgeType.equals("Ca")) {
			atomNames = new String[1];
			atomNames[0] = "CA";
		}
		// TODO add other cases for edgeType "Cb", "ALL", etc.
		
		AtomContactSet atomContacts = StructureTools.getAtomsInContact(chain, atomNames, cutoff);
		
		GroupContactSet groupContacts = new GroupContactSet(atomContacts);
		
		RIGraph graph = new RIGraph();
		
		graph.setCutoff(cutoff);
		graph.setContactType(edgeType);
		graph.setPdbCode(chain.getStructure().getPDBCode());
		graph.setChainCode(chain.getId());
		graph.setPdbChainCode(chain.getName());
		graph.setModel(modelSerial-1);
		// TODO check if we need seqres or atom sequence
		graph.setSequence(chain.getSeqResSequence());
		// TODO secondary structure
		//graph.setSecondaryStructure(this.secondaryStructure);
		
		List<Group> seqresgroups = chain.getSeqResGroups();
		
		for (Group group : seqresgroups) {
			graph.addVertex(new RIGNode(seqresgroups.indexOf(group)+1, group.getPDBName()));
		}
		
		for (GroupContact contact : groupContacts) {
			Group igroup = contact.getPair().getFirst();
			Group jgroup = contact.getPair().getSecond();
			graph.addEdgeIJ(seqresgroups.indexOf(igroup)+1, seqresgroups.indexOf(jgroup)+1);
		}
		
		return graph;
	}
	
	

	
	
	
	public void init(Sequence sequence, RIGraph contacts, SecondaryStructure ss, Chain coordinates, MySQLConnection conn) {
		//this.sequence = sequence;
		//this.contacts = contacts;
		
	}
	

	
	public static int getResSerial(ResidueNumber resNum, Chain pdb) {
		try {
			Group g = pdb.getGroupByPDB(resNum);
			return pdb.getEntityInfo().getAlignedResIndex(g, pdb);
		} catch (StructureException e) {
			System.err.println("Warning: no group found for "+resNum);
			return -1;
		}
	}
	
	public static SecondaryStructure convertSecondStruc(Structure fullpdb, Chain pdb){
		
		
		SecStrucCalc ssp = new SecStrucCalc();
		try{
			ssp.calculate(fullpdb, true);
		}
		catch(StructureException e){
			System.err.println("Warning: Cannot calculate and assign secondarystructure ");
		}
		
		

		
		List<org.biojava.nbio.structure.secstruc.SecStrucElement> sse = SecStrucTools.getSecStrucElements(fullpdb);
		//SecondaryStructure secondaryStructure = new SecondaryStructure(pdb.getSequence().getSeq());
		String seq = pdb.getSeqResSequence();
		SecondaryStructure owlss = new SecondaryStructure(seq);
		int Ecount = 0;
		int Hcount = 0;
		int Lcount = 0;
		for(org.biojava.nbio.structure.secstruc.SecStrucElement e: sse){
			if (e.getRange().getStart().getChainName().equals(pdb.getName()) &&
					e.getRange().getEnd().getChainName().equals(pdb.getName())) {
				//interval 
				int start = getResSerial(e.getRange().getStart(), pdb);
				int end = getResSerial(e.getRange().getEnd(), pdb);
				//Interval range = Interval(start, end);

				//type and ID
				char type;
				String id;
				if(e.getType().toString().equals("E")){
					type = 'S';
					String typeStr = "S";
					String intStr = Integer.toString(Ecount);
					id = typeStr.concat(intStr);
					Ecount++;
				}
				else if(e.getType().toString().equals("H") ||
						e.getType().toString().equals("G") ||
						e.getType().toString().equals("I")){
					type = 'H';
					String typeStr = "H";
					String intStr = Integer.toString(Hcount);
					id = typeStr.concat(intStr);
					Hcount++;
				}
				else{
					type = 'L';
					String typeStr = "L";
					String intStr = Integer.toString(Lcount);
					id = typeStr.concat(intStr);
					Lcount++;
				}

				SecStrucElement element = new SecStrucElement(type, start, end, id);
				owlss.add(element);


			}
		}
		return owlss;
	}
	
	
	
	public static HashMap<Pair<Integer>, Double> calcDistMatrixbjv(Chain pdb){
		HashMap<Pair<Integer>,Double> distMatrixRes = new HashMap<Pair<Integer>, Double>();
		for(Group gi: pdb.getAtomGroups()){
			//getting gi atom information
			org.biojava.nbio.structure.Atom i = gi.getAtom("CA");
			Point3d icoord = i.getCoordsAsPoint3d();
			ResidueNumber resgi = gi.getResidueNumber();
			int resnumgi = getResSerial(resgi, pdb);
			
			for(Group gj: pdb.getAtomGroups()){
				//getting gj atom information
				org.biojava.nbio.structure.Atom j = gj.getAtom("CA");
				Point3d jcoord = j.getCoordsAsPoint3d();
				ResidueNumber resgj = gj.getResidueNumber();
				int resnumgj = getResSerial(resgj, pdb);
				
				//calculating the distance between 
				double distance = icoord.distance(jcoord);
				
				//storing the lower left bottom information into the HashMap
				if(resnumgj > resnumgi){
					Pair<Integer> edge = new Pair<Integer>(resnumgi,resnumgj);
					distMatrixRes.put(edge, distance);	
				}
			}
		}
		return distMatrixRes;
	}
	
	
	
	
	
	
	/*public static HashMap<Pair<Integer>, Double> calcDistMatrix(String ct){
		HashMap<Pair<Integer>,Double> distMatrixAtoms = calcAtomDistMatrix(ct);

		 // mapping atom serials to residue serials
		 // TODO: we could integrate this with the code in calcAtomDistMatrix to avoid storing two distance maps in memory
		HashMap<Pair<Integer>,Double> distMatrixRes = new HashMap<Pair<Integer>, Double>();
		for (Pair<Integer> cont: distMatrixAtoms.keySet()){
			int i_resser = getResSerFromAtomSer(cont.getFirst());
			int j_resser = getResSerFromAtomSer(cont.getSecond());
			Pair<Integer> edge = new Pair<Integer>(i_resser,j_resser);
			if (distMatrixRes.containsKey(edge)) {
				distMatrixRes.put(edge, Math.min(distMatrixRes.get(edge), distMatrixAtoms.get(cont)));
			} else {
				distMatrixRes.put(edge, distMatrixAtoms.get(cont));
			}
		}

		return distMatrixRes;
	}*/

	
	
	
	
	
	
	
	
	
}
