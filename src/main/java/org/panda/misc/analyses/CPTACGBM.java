package org.panda.misc.analyses;

import org.panda.misc.causalpath.CausalPathSubnetwork;
import org.panda.utility.StreamDirection;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class CPTACGBM
{
	public static void main(String[] args) throws IOException
	{
//		generateEnrichedSubsets();
//		generateCOSMICUpstream();
//		takeASubgraph();
		generateSubgraphsForCorrelation();
	}

	public static void generateEnrichedSubsets() throws IOException
	{
		CausalPathSubnetwork.generateNeighborhoodSubgraphsForSignificantsRecursively("/home/ozgun/Analyses/CPTAC-GBM", 0.1);
	}

	public static void generateSubgraphsForCorrelation() throws IOException
	{
		CausalPathSubnetwork.writeSignifNeighForCorrBased("/home/ozgun/Analyses/CPTAC-GBM/tumors_correlation", StreamDirection.DOWNSTREAM, 0.1);
		CausalPathSubnetwork.writeGOINeighForCorrBased("/home/ozgun/Analyses/CPTAC-GBM/tumors_correlation", COSMIC, StreamDirection.UPSTREAM);
	}

	public static void generateCOSMICUpstream() throws IOException
	{
		CausalPathSubnetwork.writeGOINeighForCompBasedRecursively("/home/ozgun/Analyses/CPTAC-GBM", COSMIC, StreamDirection.UPSTREAM);
	}

	public static void takeASubgraph() throws IOException
	{
		CausalPathSubnetwork.writeGOINeighForCompBased("/home/ozgun/Analyses/CPTAC-GBM/mesenchymal_vs_others", new HashSet<>(Arrays.asList("SERPINE1", "SERPINA1")), StreamDirection.BOTHSTREAM, "causative-SERPINE1");
	}

	public static final Set<String> COSMIC = new HashSet<>(Arrays.asList(("DAXX\n" +
		"ERBB2\n" +
		"GOPC\n" +
		"HIF1A\n" +
		"IDH1\n" +
		"IDH2\n" +
		"LZTR1\n" +
		"MDM4\n" +
		"PDGFRA\n" +
		"PIK3CA\n" +
		"PIK3R1\n" +
		"PTPRD\n" +
		"ROS1\n" +
		"SALL4\n" +
		"STAG2\n" +
		"TERT\n" +
		"AKT3\n" +
		"ATRX\n" +
		"ZNF429\n").split("\n")));
}
