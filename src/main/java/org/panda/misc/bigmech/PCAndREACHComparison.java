package org.panda.misc.bigmech;

import org.panda.resource.network.SignedPC;
import org.panda.resource.network.SignedREACH;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.SiteSpecificGraph;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class PCAndREACHComparison
{
	public static void main(String[] args) throws IOException
	{
		printVenn();
		System.out.println();
		printVennWithSites();
		System.out.println();
		printCPResultVenn();
	}

	public static void printVenn()
	{
		DirectedGraph pP = SignedPC.get().getGraph(SignedType.PHOSPHORYLATES);
		DirectedGraph dP = SignedPC.get().getGraph(SignedType.DEPHOSPHORYLATES);
		DirectedGraph pR = SignedREACH.get().getGraph(SignedType.PHOSPHORYLATES);
		DirectedGraph dR = SignedREACH.get().getGraph(SignedType.DEPHOSPHORYLATES);

		Set<String> setPC = new HashSet<>();
		collectSignatures(pP, setPC);
		collectSignatures(dP, setPC);

		Set<String> setREACH = new HashSet<>();
		collectSignatures(pR, setREACH);
		collectSignatures(dR, setREACH);

		CollectionUtil.printNameMapping("PC", "REACH");
		CollectionUtil.printVennCounts(setPC, setREACH);
	}

	public static void printVennWithSites()
	{
		SiteSpecificGraph pP = (SiteSpecificGraph) SignedPC.get().getGraph(SignedType.PHOSPHORYLATES);
		SiteSpecificGraph dP = (SiteSpecificGraph) SignedPC.get().getGraph(SignedType.DEPHOSPHORYLATES);
		SiteSpecificGraph pR = (SiteSpecificGraph) SignedREACH.get().getGraph(SignedType.PHOSPHORYLATES);
		SiteSpecificGraph dR = (SiteSpecificGraph) SignedREACH.get().getGraph(SignedType.DEPHOSPHORYLATES);

		Set<String> setPC = new HashSet<>();
		collectSignaturesWithSite(pP, setPC);
		collectSignaturesWithSite(dP, setPC);

		Set<String> setREACH = new HashSet<>();
		collectSignaturesWithSite(pR, setREACH);
		collectSignaturesWithSite(dR, setREACH);

		CollectionUtil.printNameMapping("PC", "REACH");
		CollectionUtil.printVennCounts(setPC, setREACH);
	}

	private static void collectSignatures(DirectedGraph graph, Set<String> set)
	{
		graph.getOneSideSymbols(true).forEach(s ->
		{
			for (String t : graph.getDownstream(s))
			{
				set.add(s + "\t" + graph.getEdgeType() + "\t" + t);
			}
		});
	}

	private static void collectSignaturesWithSite(SiteSpecificGraph graph, Set<String> set)
	{
		graph.getOneSideSymbols(true).forEach(s ->
		{
			for (String t : graph.getDownstream(s))
			{
				Set<String> sites = graph.getSites(s, t);
				if (sites != null)
				{
					for (String site : sites)
					{
						set.add(s + "\t" + graph.getEdgeType() + "\t" + t + "\t" + site);
					}
				}
			}
		});
	}

	private static void printCPResultVenn() throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/BigMech/Diff-view/";
		Set<String> setPC = loadCPResult(dir + "pc.sif");
		Set<String> setREACH = loadCPResult(dir + "reach.sif");
		CollectionUtil.printNameMapping("PC", "REACH");
		CollectionUtil.printVennCounts(setPC, setREACH);
	}

	private static Set<String> loadCPResult(String file) throws IOException
	{
		return Files.lines(Paths.get(file)).map(l -> l.split("\t")).filter(t -> t.length > 2)
			.map(t -> ArrayUtil.getString("\t", t[0], t[1], t[2])).collect(Collectors.toSet());
	}
}
