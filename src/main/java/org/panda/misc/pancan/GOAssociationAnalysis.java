package org.panda.misc.pancan;

import org.panda.resource.GO;
import org.panda.resource.HGNC;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.statistics.FishersExactTest;
import org.panda.utility.statistics.Overlap;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class GOAssociationAnalysis
{
	public static void main(String[] args) throws IOException
	{
		printResultGeneStatusForHighlighting();
//		printNegativeControlGeneStatusForHighlighting();
//		printEnrichmentMovingWindow();
	}

	static void printResultGeneStatusForHighlighting() throws IOException
	{
		int size = 100;
		int cols = 5;
		List<String> genes = readGeneList(size);

		Set<String> ids = readSelectedEnrichedGOIDs();
		Set<String> goTagged = new HashSet<>(GO.get().getGenes(ids));
		System.out.println("goTagged.size() = " + goTagged.size());

		Set<String> canTagged = new HashSet<>(PanCanResultLoader.readCancerGenes());
		System.out.println("canTagged.size() = " + canTagged.size());

		Set<String> pathwayEnriched = Files.lines(Paths.get("/home/babur/Documents/PanCan/tissue-unnormalized-results/pathway-enrichment.txt"))
			.skip(10).map(l -> l.split("\t")).filter(t -> Double.parseDouble(t[4]) < 0.01).map(t -> t[8].split(", "))
			.flatMap(Arrays::stream).collect(Collectors.toSet());

		System.out.println("pathwayEnriched.size() = " + pathwayEnriched.size());

		System.out.println("All genes in GO = " + GO.get().getAllGenes().size());


		List<String> tag = new ArrayList<>(size);
		for (String gene : genes)
		{
			tag.add(canTagged.contains(gene) ? "X" : goTagged.contains(gene) ? pathwayEnriched.contains(gene) ? "EB" : "EG" : pathwayEnriched.contains(gene) ? "EP" : "");
		}
		System.out.println("X = " + CollectionUtil.countValues(tag, "X"));
		System.out.println("EG = " + CollectionUtil.countValues(tag, "EG"));
		System.out.println("EP = " + CollectionUtil.countValues(tag, "EP"));
		System.out.println("EB = " + CollectionUtil.countValues(tag, "EB"));
		System.out.println("no tag = " + CollectionUtil.countValues(tag, ""));
		printInColumns(genes, size/cols);
		printInColumns(tag, size/cols);
	}

	static void printNegativeControlGeneStatusForHighlighting() throws IOException
	{
		int size = 500;
		int cols = 10;
		List<String> genes = readControlList(size);

		Set<String> ids = readSelectedEnrichedControlGOIDs();
		Set<String> goTagged = new HashSet<>(GO.get().getGenes(ids));
		System.out.println("goTagged.size() = " + goTagged.size());

		Set<String> canTagged = new HashSet<>(PanCanResultLoader.readCancerGenes());
		System.out.println("canTagged.size() = " + canTagged.size());

		Set<String> pathwayEnriched = Files.lines(Paths.get("/home/babur/Documents/PanCan/tissue-unnormalized-results/pathway-enrichment-control.txt"))
			.skip(10).map(l -> l.split("\t")).filter(t -> Double.parseDouble(t[4]) < 0.01).map(t -> t[8].split(", "))
			.flatMap(Arrays::stream).collect(Collectors.toSet());

		System.out.println("pathwayEnriched.size() = " + pathwayEnriched.size());

		System.out.println("All genes in GO = " + GO.get().getAllGenes().size());


		List<String> tag = new ArrayList<>(size);
		for (String gene : genes)
		{
			tag.add(canTagged.contains(gene) ? "X" : goTagged.contains(gene) ? pathwayEnriched.contains(gene) ? "EB" : "EG" : pathwayEnriched.contains(gene) ? "EP" : "");
		}
		System.out.println("X = " + CollectionUtil.countValues(tag, "X"));
		System.out.println("EG = " + CollectionUtil.countValues(tag, "EG"));
		System.out.println("EP = " + CollectionUtil.countValues(tag, "EP"));
		System.out.println("EB = " + CollectionUtil.countValues(tag, "EB"));
		System.out.println("no tag = " + CollectionUtil.countValues(tag, ""));
		printInColumns(genes, size/cols);
		printInColumns(tag, size/cols);
	}

	static Set<String> readInterestKeywords() throws IOException
	{
		return Files.lines(Paths.get("/home/babur/Documents/PanCan/GO-term-interest.txt"))
			.filter(l -> !l.startsWith("#")).collect(Collectors.toSet());
	}

	static Set<String> readSelectedEnrichedGOIDs() throws IOException
	{
		return Files.lines(Paths.get("/home/babur/Documents/PanCan/selected-enriched-GO-terms.txt"))
			.filter(l -> !l.startsWith("#")).map(l -> l.split("\t")[0]).collect(Collectors.toSet());
	}

	static Set<String> readSelectedEnrichedControlGOIDs() throws IOException
	{
		return Files.lines(Paths.get("/home/babur/Documents/PanCan/selected-enriched-GO-terms-control.txt"))
			.filter(l -> !l.startsWith("#")).map(l -> l.split("\t")[0]).collect(Collectors.toSet());
	}

	static List<String> readGeneList(int limit) throws IOException
	{
		return Files.lines(Paths.get("/home/babur/Documents/PanCan/pancan.txt")).skip(1).limit(limit)
//		return Files.lines(Paths.get("/home/babur/Documents/PanCan/tissue-unnormalized-results/pancan.txt")).skip(1).limit(limit)
			.map(l -> l.split("\t")[0]).collect(Collectors.toList());
	}

	static List<String> readControlList(int limit) throws IOException
	{
		return Files.lines(Paths.get("/home/babur/Documents/PanCan/tissue-unnormalized-results/sorted-to-freq.txt")).skip(1).limit(limit)
			.map(l -> l.split("\t")[0]).collect(Collectors.toList());
	}

	static <T> void printInColumns(List<T> list, int rows)
	{
		for (int i = 0; i < rows; i++)
		{
			for (int j = i; j < list.size(); j+=rows)
			{
				System.out.print(list.get(j) + "\t");
			}
			System.out.println();
		}
	}

	static void printEnrichmentMovingWindow() throws IOException
	{
		List<String> genes = readGeneList(Integer.MAX_VALUE);
		List<String> freqSorted = getSortedToMutationFreq(genes);

//		Set<String> allGOGenes = new HashSet<>(GO.get().getAllGenes());
//		Set<String> assocGenes = GO.get().getGenesContainingKeywordInTermNames(readInterestKeywords());

		// Use cancer genes instead of GO terms
		Set<String> allGOGenes = new HashSet<>(genes);
		Set<String> assocGenes = PanCanResultLoader.readCancerGenes();

		int window = 100;

		allGOGenes.retainAll(genes);
		assocGenes.retainAll(genes);

		double expectedRatio = assocGenes.size() / (double) allGOGenes.size();
		System.out.println(0 + "\t" + expectedRatio + "\texpected ratio");
		System.out.println(7000 + "\t" + expectedRatio);

		System.out.print("\nFrom\tTo\tGenes\tGenes in GO\tAssociated\tRatio\tp-val\tCtrl ratio\tCtrl p-val");
		for (int i = 0; i < genes.size() - window; i+=10)
		{
			Set<String> sub = new HashSet<>(genes.subList(i, i + window));

			int readCnt = sub.size();
			sub.retainAll(allGOGenes);
			int readInGO = sub.size();
			sub.retainAll(assocGenes);
			int ass = sub.size();

			System.out.print("\n" + i + "\t" + (i + window) + "\t" + readCnt + "\t" + readInGO + "\t" + ass + "\t" +
				(ass / (double) readInGO) + "\t" + FishersExactTest.calcEnrichmentPval(allGOGenes.size(), assocGenes.size(), readInGO, ass));

			sub = new HashSet<>(freqSorted.subList(i, i + window));

			sub.retainAll(allGOGenes);
			readInGO = sub.size();
			sub.retainAll(assocGenes);
			ass = sub.size();

			System.out.print("\t" + (ass / (double) readInGO) + "\t" +
				FishersExactTest.calcEnrichmentPval(allGOGenes.size(), assocGenes.size(), readInGO, ass));
		}
	}

	public static List<String> getSortedToMutationFreq(List<String> original) throws IOException
	{
		List<String> genes = new ArrayList<>(original);
		Map<String, Integer> altCnt = new HashMap<>();
		Files.lines(Paths.get("/home/babur/Documents/mutex/TCGA/PanCan/1/1/DataMatrix.txt")).skip(1).forEach(l ->
		{
			String gene = l.substring(0, l.indexOf("\t"));
			if (genes.contains(gene))
			{
				String[] vals = l.substring(l.indexOf("\t") + 1).split("\t");
				int count = ArrayUtil.countValues(vals, "1");
				altCnt.put(gene, count);
			}
		});

		Collections.sort(genes, (o1, o2) -> altCnt.get(o2).compareTo(altCnt.get(o1)));
		return genes;
	}
}
