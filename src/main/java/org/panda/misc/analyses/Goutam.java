package org.panda.misc.analyses;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.misc.PathwayEnrichmentSIFGenerator;
import org.panda.resource.PCPathway;
import org.panda.utility.CollectionUtil;
import org.panda.utility.Kronometre;
import org.panda.utility.statistics.FDR;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.*;

/**
 * Created by babur on 2/29/16.
 */
public class Goutam
{
	static String DIR = "/home/babur/Documents/Analyses/Goutam/";

	public static void main(String[] args) throws IOException
	{
		Kronometre k = new Kronometre();

		Set<String> allGenes = getUnion(readAllGenes(0));
		Set<String> upGenes = getUnion(readAllGenes(1));
		Set<String> dwGenes = getUnion(readAllGenes(-1));
		Set<String> undecidedGenes = CollectionUtil.getIntersection(upGenes, dwGenes);

		Set<String> pathwayIDs = doEnrichment(1);
		pathwayIDs.addAll(doEnrichment(-1));

		PathwayEnrichmentSIFGenerator sg = new PathwayEnrichmentSIFGenerator();
		sg.setGenes(allGenes);
		sg.setOwlFilename("/home/babur/Documents/PC/PathwayCommons.8.Detailed.BIOPAX.owl");
		sg.setBlacklistFilename("/home/babur/Documents/PC/blacklist_v8.txt");
		sg.setPathwayIDs(pathwayIDs);
		String[][] groupTerms = {
			new String[]{"collagen", "extracellular"},
			new String[]{"integrin"},
			new String[]{"AP-1", "AP1"},
//			new String[]{"TGF-beta", "TGF_beta"},
//			new String[]{"syndecan"},
//			new String[]{"fibrin"},
//			new String[]{"RA biosynthesis", "retinoid", "visual"},
//			new String[]{"C-MYC"},
//			new String[]{"HIF-1"},
//			new String[]{"Scavenging"},
//			new String[]{"VEGF"},
//			new String[]{"FGF"},
//			new String[]{"miRNA"},
//			new String[]{"pyrimidine"},
//			new String[]{"Creatine"}
		};

		sg.setGroupTerms(groupTerms);
		sg.setTypes(SIFEnum.CONTROLS_STATE_CHANGE_OF, SIFEnum.CONTROLS_EXPRESSION_OF);

		Color upColor = new Color(180, 255, 200);
		Color dwColor = new Color(255, 180, 200);
		Color mixColor = new Color(255, 255, 220);
		allGenes.forEach(g -> sg.addNodeColor(g,
			undecidedGenes.contains(g) ? mixColor : upGenes.contains(g) ? upColor : dwColor));

		sg.write(DIR + "pathway2");
		k.print();
	}

	private static Set<String> doEnrichment(int sign) throws IOException
	{
		System.out.println("sign = " + sign);
		List<Set<String>> sets = readAllGenes(sign);
		Set<String> genes = getIntersection(sets);
		showOverlaps(sets);
		return findEnrichedPathways(genes, sign);
	}

	private static void showOverlaps(List<Set<String>> bag) throws IOException
	{
		CollectionUtil.printNameMapping("Caso Cl 2", "MDV Cl 1", "Neo Sh#5", "Neo Sh# 3");
		CollectionUtil.printVennSets(bag.toArray(new Collection[bag.size()]));
	}

	private static Set<String> findEnrichedPathways(Set<String> genes, Integer sign) throws IOException
	{
		String filename = DIR + "results" + (sign == 0 ? "" : sign > 0 ? "-upregulated" : "-downregulated") + ".txt";

		int minMember = 3;
		int maxMember = 300;
		PCPathway.get().writeEnrichmentResults(genes, minMember, maxMember, filename);
		Map<String, Double>[] pvals = PCPathway.get().getEnrichmentPvals(genes, null, minMember, maxMember);
		Map<String, Double> qvals = FDR.getQVals(pvals[0], pvals[1]);
		OptionalDouble thr = pvals[0].keySet().stream().filter(id -> qvals.get(id) < 0.1).mapToDouble(pvals[0]::get).max();
		if (!thr.isPresent()) return null;

		return pvals[0].keySet().stream().filter(id -> pvals[0].get(id) <= thr.getAsDouble()).collect(toSet());
	}

	private static Set<String> getIntersection(List<Set<String>> bag)
	{
		Set<String> intersection = new HashSet<>();

		boolean first = true;

		for (Set<String> set : bag)
		{
			if (first)
			{
				intersection.addAll(set);
				first = false;
			}
			else
			{
				intersection.retainAll(set);
			}
		}
		return intersection;
	}

	private static Set<String> getUnion(List<Set<String>> bag)
	{
		return bag.stream().flatMap(Collection::stream).collect(Collectors.toSet());
	}

	private static List<Set<String>> readAllGenes(Integer sign) throws FileNotFoundException
	{
		List<Set<String>> bag = new ArrayList<>();
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/Caso Cl 2 vs Contol copy.csv", sign));
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/MDV Cl 1 vs Control copy.csv", sign));
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/Neo Sh#5  vs Control copy.csv", sign));
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/Neo Sh# 3 vs Control copy.csv", sign));
		return bag;
	}

	private static Set<String> readGenes(String file, int sign) throws FileNotFoundException
	{
		Set<String> set = new HashSet<>();
		Scanner sc = new Scanner(new File(file));
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			String gene = token[2];
			if (gene.isEmpty()) continue;
			String foldChange = token[5];

			if (sign == 0 || (sign < 0 && foldChange.startsWith("-")) || sign > 0 && !foldChange.startsWith("-"))
				set.add(gene);
		}
		return set;
	}


}
