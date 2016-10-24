package cmview.datasources;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;

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
import owl.core.sequence.alignment.MultipleSequenceAlignment;
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
	
	
	
	
	public static boolean containstdAaResidueBioj(int resSerial, Chain pdb){
		if(pdb.getSeqResGroup(resSerial-1)!=null && pdb.getSeqResGroup(resSerial-1).isAminoAcid()){
			Group g = pdb.getSeqResGroup(resSerial-1);
			if (pdb.getAtomGroups().contains(g)) {
				return true;
			} 
		}
		
		return false;
	}
	
	public static TreeMap<Integer, Character> getUnobservedResiduesBioj(Chain c){
		TreeMap<Integer, Character> tmap = new TreeMap<Integer, Character>();
		for(int i = 1; i < c.getSeqResLength(); i++){
			if(!containstdAaResidueBioj(i,c)){
				//how to get the residue?
				Character value = StructureTools.get1LetterCodeAmino(c.getAtomGroup(i).getPDBName());
				tmap.put(i, value);
			}
		}
		return tmap;
	}
	
	
	
	public static int aAResidueCountBioj(Chain pdb){
		int count = 0;
		for(int i = 0; i < pdb.getAtomLength(); i++){
			if(pdb.getAtomGroup(i).isAminoAcid()){
				count++;
			}
		}
		return count;
	}
	

	//contacttype = Ca
	public static HashMap<Pair<Integer>,Double> getDiffDistMapbiojchain(String contactType1, Chain pdb1,
			Chain pdb2, String contactType2, MultipleSequenceAlignment ali, String name1, String name2){
		HashMap<Pair<Integer>,Double> otherDistMatrix = calcDistMatrixbjv(pdb2);
		HashMap<Pair<Integer>,Double> thisDistMatrix = calcDistMatrixbjv(pdb1);
		//get sequence length = getSeqResLength?
		//System.out.println("pdb1.length: " + pdb1.getSeqResLength() + " pdb2.length: " + pdb2.getSeqResLength());
		HashMap<Pair<Integer>,Double> alignedDistMatrix = new HashMap<Pair<Integer>, Double>(Math.min(pdb1.getSeqResLength(), pdb2.getSeqResLength()));
		
		int i1,i2,j1,j2;
		TreeSet<Integer> unobserved1 = new TreeSet<Integer>();
		TreeSet<Integer> unobserved2 = new TreeSet<Integer>();
		// detect all unobserved residues
		
		
		/*
		 * containsStdAaResidue:
		 * Tells whether residue of given residue number is a (observed) standard amino 
		 * acid residue in this PdbChain instance.  
		 * 
		 */
		for(int i = 1; i <= ali.getAlignmentLength(); ++i) {
			i1 = ali.al2seq(name1, i);
			i2 = ali.al2seq(name2, i);
			//containsStdAaResidue(i1)
			if( i1 != -1 &&	!containstdAaResidueBioj(i1,pdb1)) {
				unobserved1.add(i1);
			}
			
			if( i2 != -1 && !containstdAaResidueBioj(i2,pdb2)) {
				unobserved2.add(i2);
			}
		}
		
		
		
		// strategy: we always have to look through the alignment to say 
		// whether a difference distance can be assigned to a pair of 
		// corresponding residues. To put it differently, for any two 
		// alignment columns one always has to ensure that both columns 
		// only contain observed residues (no gaps!), otherwise the one 
		// cannot obtain a distance in at least one structure as a gap 
		// indicates "no coordinates available".  

		for(int i = 1; i <= ali.getAlignmentLength()-1; ++i) {

			i1 = ali.al2seq(name1, i);
			i2 = ali.al2seq(name2, i);

			// alignment columns must not contain gap characters and both 
			// residues in the current column have to be observed!
			if( i1 == -1 || i2 == -1 || unobserved1.contains(i1) || unobserved2.contains(i2) ) {
				continue;
			}

			for(int j = i + 1; j <= ali.getAlignmentLength(); ++j) {

				j1 = ali.al2seq(name1, j);
				j2 = ali.al2seq(name2, j);

				if( j1 == -1 || j2 == -1 || unobserved1.contains(j1) || unobserved2.contains(j2) ) {
					continue;
				}

				// make the edges
				Pair<Integer> e1 = new Pair<Integer>(i1,j1);
				Pair<Integer> e2 = new Pair<Integer>(i2,j2);

				if(thisDistMatrix.get(e1) == null || otherDistMatrix.get(e2) == null){
					System.err.println("Warning: edge value equals null, e1: " + e1 + " e2: " + e2);
					continue;
				}
				alignedDistMatrix.put(new Pair<Integer>(i,j),Math.abs(thisDistMatrix.get(e1)-otherDistMatrix.get(e2)));
			}
		}
		return alignedDistMatrix;
	}
	
	
	
	

	public static void writeToPDBFileBioj(File outFile, Chain c) throws FileNotFoundException{
		PrintStream out = new PrintStream(new FileOutputStream(outFile));
		String pdboutput = c.toPDB();
		out.println(pdboutput);
		
		out.println("END");
		out.close();
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}
