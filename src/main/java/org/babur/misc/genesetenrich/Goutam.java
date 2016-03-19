package org.babur.misc.genesetenrich;

import org.cbio.causality.network.PCPathway;
import org.cbio.causality.util.CollectionUtil;

import java.io.*;
import java.util.*;

/**
 * Created by babur on 2/29/16.
 */
public class Goutam
{
	public static void main(String[] args) throws IOException
	{
		Integer sign = 1;
		System.out.println("sign = " + sign);
		List<Set<String>> sets = readAllGenes(sign);
		Set<String> genes = getIntersection(sets);
		showOverlaps(sets);
		findEnrichedPathways(genes, sign);
	}

	private static void showOverlaps(List<Set<String>> bag) throws IOException
	{
		CollectionUtil.printNameMapping("Caso Cl 2", "MDV Cl 1", "Neo Sh#5", "Neo Sh# 3");
		CollectionUtil.printVennSets(bag.toArray(new Collection[bag.size()]));
	}

	private static void findEnrichedPathways(Set<String> genes, Integer sign) throws IOException
	{
		String filename = "/home/babur/Downloads/AnalysisList/results" + (sign == null ? "" : sign > 0 ? "-upregulated"
			: "-downregulated") + ".txt";

		PCPathway.writeEnrichmentResults(genes, 3, filename);
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

	private static List<Set<String>> readAllGenes(Integer sign) throws FileNotFoundException
	{
		List<Set<String>> bag = new ArrayList<>();
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/Caso Cl 2 vs Contol copy.csv", sign));
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/MDV Cl 1 vs Control copy.csv", sign));
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/Neo Sh#5  vs Control copy.csv", sign));
		bag.add(readGenes("/home/babur/Downloads/AnalysisList/Neo Sh# 3 vs Control copy.csv", sign));
		return bag;
	}

	private static Set<String> readGenes(String file, Integer sign) throws FileNotFoundException
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

			if (sign == null || (sign < 0 && foldChange.startsWith("-")) || sign > 0 && !foldChange.startsWith("-"))
				set.add(gene);
		}
		return set;
	}


}
