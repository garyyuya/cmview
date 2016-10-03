package cmview.datasources;

import java.util.List;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.contact.GroupContact;
import org.biojava.nbio.structure.contact.GroupContactSet;

import owl.core.sequence.Sequence;
import owl.core.structure.features.SecondaryStructure;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
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
}
