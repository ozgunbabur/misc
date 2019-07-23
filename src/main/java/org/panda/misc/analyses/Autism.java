package org.panda.misc.analyses;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.panda.resource.autismdatasets.DenovoDB;
import org.panda.resource.ReactomePathway;
import org.panda.resource.autismdatasets.SFARI;
import org.panda.utility.*;
import org.panda.utility.statistics.*;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class Autism
{
	public static void main(String[] args) throws IOException
	{
//		checkEnrichmentForSFARIGenes();
//		printBiasToAutism();

//		compareUniformityVErsusFreqROCs();

//		getPvaluesForBiasToYuenDataset();
//		Map<String, Double>[] p = getPvaluesForBiasToAutism();
//		generateAlterationMatrixFromDenovoDB(p[0].keySet().stream().sorted(Comparator.comparing(p[0]::get)).limit(100).collect(Collectors.toSet()));

//		testMutexnessOfSFARI();
//		testMutexnessOfReactomePathways();

//		printOverlapsOfMutations();
//		selectFromMutexResults();
//		printSelectPathways();
//		printSFARIRanksOfGenes();
//		findSignificantMembersOfGroups();

//		printBiasOfAMutexGroup();
		printBiasOfMutexGroups();
//		printBiasOfReactomePathways();
//		printDatasetBiasOfReactomePathways();

//		compareTurnerAndAn();
	}


	/**
	 * Reads autism mutations from denovo-db, prints ranked list of genes according to their number of mutations.
	 */
	private static List<String> rankToNumberOfMutationsInDenovoDB(Map<String, Set<String>> hitMap)
	{
		return hitMap.keySet().stream()
//			.filter(g -> hitMap.get(g).size() >= 5)
			.sorted((o1, o2) -> Integer.compare(hitMap.get(o2).size(), hitMap.get(o1).size()))
//			.peek(g -> System.out.println(g + "\t" + hitMap.get(g).size()))
			.collect(Collectors.toList());
	}


	private static Map<String, Set<String>> readHitMapForAutismFromDenovoDB()
	{
		DenovoDB.DataFilter filter = e -> DenovoDB.SELECT_MUT.select(e) && DenovoDB.SELECT_AUTISM_STUDY.select(e) &&
			e.primaryPhenotype.equals("autism");

		return EQUALIZE_MUTS ? readHitMap(filter, MUT_CNT) : readHitMap(filter);
	}

	private static Map<String, Set<String>> readHitMapForGONLFromDenovoDB()
	{
		DenovoDB.DataFilter filter = e -> DenovoDB.SELECT_MUT.select(e) && e.primaryPhenotype.equals("control") &&
			e.studyName.equals("GONL");

		return EQUALIZE_MUTS ? readHitMap(filter, MUT_CNT) : readHitMap(filter);
	}

	private static Map<String, Set<String>> readHitMapForControlsInAutismStudiesFromDenovoDB()
	{
		DenovoDB.DataFilter filter = e -> DenovoDB.SELECT_MUT.select(e) && DenovoDB.SELECT_AUTISM_STUDY.select(e) &&
			e.primaryPhenotype.equals("control");

		return EQUALIZE_MUTS ? readHitMap(filter, MUT_CNT) : readHitMap(filter);
	}

	private static Map<String, Set<String>> readHitMapForOtherDisorderFromDenovoDB()
	{
		DenovoDB.DataFilter filter = e -> DenovoDB.SELECT_MUT.select(e) && !e.primaryPhenotype.equals("autism") &&
			!e.primaryPhenotype.equals("control") && !e.primaryPhenotype.equals("mixed") &&
			!e.primaryPhenotype.equals("developmentalDisorder");

		return EQUALIZE_MUTS ? readHitMap(filter, MUT_CNT) : readHitMap(filter);
	}

	private static Map<String, Set<String>> readHitMap(DenovoDB.DataFilter filter)
	{
//		Map<String, Long> cnt = getPerSampleMutationsHistogram(filter);
		Map<String, Set<String>> hitMap = new HashMap<>();
		DenovoDB.get().getDataStream(filter).forEach(e ->
		{
			if (!hitMap.containsKey(e.gene)) hitMap.put(e.gene, new HashSet<>());
			hitMap.get(e.gene).add(e.sampleID);
		});
//		getMutationCountInHitMap(hitMap);
		return hitMap;
	}

	private static int getMutationCountInHitMap(Map<String, Set<String>> hitMap)
	{
		int[] cc = new int[]{0};
		hitMap.values().stream().map(Set::size).forEach(c -> cc[0] += c);
		System.out.println("mut count in hitMap = " + cc[0]);
		return cc[0];
	}

	private static Map<String, Set<String>> readHitMap(DenovoDB.DataFilter filter, int limit)
	{
		Map<String, Set<String>> hitMap = readHitMap(filter);
		hitMap = shrinkHitMap(hitMap, limit);
		getMutationCountInHitMap(hitMap);
		return hitMap;
	}

	private static Map<String, Set<String>> shrinkHitMap(Map<String, Set<String>> hitMap, int limit)
	{
		List<String> keys = new ArrayList<>();
		hitMap.forEach((gene, pos) ->
		{
			for (String p : pos)
			{
				keys.add(gene + "|" + p);
			}
		});
		Collections.shuffle(keys);
		System.out.println("keys.size() = " + keys.size());

		Map<String, Set<String>> map = new HashMap<>();
		keys.subList(0, Math.min(limit, keys.size())).forEach(k ->
		{
			String[] t = k.split("\\|");
			if (!map.containsKey(t[0])) map.put(t[0], new HashSet<>());
			map.get(t[0]).add(t[1]);
		});
		return map;
	}

	private static Map<String, Long> getPerSampleMutationsHistogram(DenovoDB.DataFilter filter)
	{
		Map<String, Long> map = DenovoDB.get().getDataStream(filter).collect(
			Collectors.groupingBy(e -> e.studyName + e.sampleID, Collectors.counting()));

		Histogram h = new Histogram(10);
		h.setBorderAtZero(true);
		map.forEach((entry, aLong) -> h.count(aLong));
		h.printDensity();

		return map;
	}

	private static final boolean EQUALIZE_MUTS = false;
	private static final int MUT_CNT = 4923;

	private static void checkEnrichmentForSFARIGenes() throws IOException
	{
//		List<String> ranked = readMutexResultsAsRankedList("/home/ozgun/Documents/Grants/Mental/mutex/on-network/ranked-groups.txt");

//		List<String> ranked = Files.lines(Paths.get(
////			"/home/ozgun/Documents/Grants/Mental/mutation-heterogenetity-ranked.txt"))
//			"/home/ozgun/Documents/Grants/Mental/number-of-mutations-ranked.txt"))
//			"/home/ozgun/Documents/Grants/Mental/number-of-mutations-in-everything-ranked.txt"))
////			"/home/ozgun/Documents/Grants/Mental/density-of-mutations-ranked.txt"))
//			.map(l -> l.split("\t")[0]).collect(Collectors.toList());

		List<String> rankedForAutism = rankToNumberOfMutationsInDenovoDB(readHitMapForAutismFromDenovoDB());
		List<String> rankedForControls = rankToNumberOfMutationsInDenovoDB(readHitMapForControlsInAutismStudiesFromDenovoDB());

		System.out.println("\nRank\tAutism AUC     \tControl AUC");
		for (int i = 1; i <= 5; i++)
		{
			Set<String> sfariGenes = new HashSet<>(SFARI.get().getGenesWithMaxScore(i));
			System.out.println(i + "\t" + ROC.getAUC(rankedForAutism, sfariGenes) + "\t" + ROC.getAUC(rankedForControls, sfariGenes));

			PrintStream out = new PrintStream(new FileOutputStream("/home/ozgun/Documents/Grants/Mental/ROC" + i + ".txt"));
			ROC.printPlotForGoogleSheets(sfariGenes, out, rankedForAutism, rankedForControls);
		}
	}

	private static List<String> printBiasToAutism()
	{
		Map<String, Double>[] pl = getPvaluesForBiasToAutism();

		List<String> select = FDR.select(pl[0], pl[1], 0.2);
		System.out.println("select.size() = " + select.size());

		List<String> ranked = pl[0].keySet().stream().sorted(Comparator.comparing(pl[0]::get)).collect(Collectors.toList());
		Set<String> sfariGenes = new HashSet<>(SFARI.get().getGenesWithMaxScore(2));
		System.out.println("AUC = " + ROC.getAUC(ranked, sfariGenes));

		select.forEach(g -> System.out.println((SFARI.get().isAutismGene(g) ? "+" : "-") + "\t" + g + "\t" + pl[0].get(g)));

		return ranked;
	}

	private static Map<String, Double>[] getPvaluesForBiasToAutism()
	{
		Map<String, Set<String>> aHit = readHitMapForAutismFromDenovoDB();
		Map<String, Set<String>> cHit = readHitMapForControlsInAutismStudiesFromDenovoDB();

		keepPairedSamplesRemoveOverlaps(aHit, cHit);

		int aTotal = getMutationCountInHitMap(aHit);
		int cTotal = getMutationCountInHitMap(cHit);

		double bias = cTotal / (double) (aTotal + cTotal);
		System.out.println("bias = " + bias);

		Map<String, Double> pvals = new HashMap<>();
		Map<String, Double> limits = new HashMap<>();

		for (String gene : aHit.keySet())
		{
			int aCnt = aHit.get(gene).size();
			int cCnt = cHit.containsKey(gene) ? cHit.get(gene).size() : 0;

			BinomialDistribution dist = new BinomialDistribution(aCnt + cCnt, bias);
			double p = dist.cumulativeProbability(cCnt);
			double l = dist.cumulativeProbability(0);

			pvals.put(gene, p);
			limits.put(gene, l);
		}

		aHit.keySet().stream()
//			.sorted((o1, o2) -> Integer.compare(aHit.get(o2).size(), aHit.get(o1).size()))
			.sorted(Comparator.comparing(pvals::get))
			.forEach(gene ->
			System.out.println(SFARI.get().getClassification(gene) + "\t" + gene + "\t" + aHit.get(gene).size() + "\t" +
				cHit.getOrDefault(gene, Collections.emptySet()).size() + "\t" + pvals.get(gene)));

		List<String> select = FDR.select(pvals, limits, 0.1);
		System.out.println("select.size() = " + select.size());


		return new Map[]{pvals, limits};
	}

	private static void keepPairedSamplesRemoveOverlaps(Map<String, Set<String>> aHit, Map<String, Set<String>> cHit)
	{
//		aHit.values().forEach(set -> set.removeIf(s -> !(s.endsWith(".p1") || s.endsWith(".s1"))));
//		cHit.values().forEach(set -> set.removeIf(s -> !(s.endsWith(".p1") || s.endsWith(".s1"))));
//		if (true) return;


		// crop to controlled cases
//		Set<String> sampleKeys = DenovoDB.get().getControlledSamples();
//		aHit.values().forEach(set -> set.retainAll(sampleKeys));
//		cHit.values().forEach(set -> set.retainAll(sampleKeys));

		// find and remove overlaps
		for (String gene : aHit.keySet())
		{
			if (cHit.containsKey(gene))
			{
				Set<String> common = getCommonFamilies(aHit, cHit, gene);
				aHit.get(gene).removeAll(common);
				cHit.get(gene).removeAll(common);
			}
		}

		// remove empty genes
		new HashSet<>(aHit.keySet()).forEach(gene ->
		{
			if (aHit.get(gene).isEmpty()) aHit.remove(gene);
		});
		new HashSet<>(cHit.keySet()).forEach(gene ->
		{
			if (cHit.get(gene).isEmpty()) cHit.remove(gene);
		});
	}

	private static Set<String> getCommonFamilies(Map<String, Set<String>> aHit, Map<String, Set<String>> cHit, String gene)
	{
		Set<String> common = new HashSet<>();

		for (String aKey : aHit.get(gene))
		{
			String cKey = DenovoDB.get().getControlOfPatient(aKey);
			if (cKey != null && cHit.containsKey(gene) && cHit.get(gene).contains(cKey))
			{
				common.add(aKey);
				common.add(cKey);
			}
		}
		return common;
	}

	private static Map<String, Double>[] getPvaluesForBiasToYuenDataset()
	{
		DenovoDB.DataFilter filter = e -> DenovoDB.SELECT_MUT.select(e) && DenovoDB.SELECT_AUTISM_STUDY.select(e) &&
			e.studyName.equals("Yuen2017") && e.primaryPhenotype.equals("autism");

		Map<String, Set<String>> cHit = readHitMap(filter);

		filter = e -> DenovoDB.SELECT_MUT.select(e) && DenovoDB.SELECT_AUTISM_STUDY.select(e) &&
			e.studyName.equals("Turner_2017") && e.primaryPhenotype.equals("autism");

		Map<String, Set<String>> aHit = readHitMap(filter);

		int aTotal = getMutationCountInHitMap(aHit);
		int cTotal = getMutationCountInHitMap(cHit);

		double bias = cTotal / (double) (aTotal + cTotal);
		System.out.println("bias = " + bias);

		Map<String, Double> pvals = new HashMap<>();
		Map<String, Double> limits = new HashMap<>();

		for (String gene : aHit.keySet())
		{
			int aCnt = aHit.get(gene).size();
			int cCnt = cHit.containsKey(gene) ? cHit.get(gene).size() : 0;

			BinomialDistribution dist = new BinomialDistribution(aCnt + cCnt, bias);
			double p = dist.cumulativeProbability(cCnt);
			double l = dist.cumulativeProbability(0);

			pvals.put(gene, p);
			limits.put(gene, l);
		}

		aHit.keySet().stream()
//			.sorted((o1, o2) -> Integer.compare(aHit.get(o2).size(), aHit.get(o1).size()))
			.sorted(Comparator.comparing(pvals::get))
			.forEach(gene ->
				System.out.println(SFARI.get().getClassification(gene) + "\t" + gene + "\t" + aHit.get(gene).size() + "\t" +
					cHit.getOrDefault(gene, Collections.emptySet()).size() + "\t" +
					(getCommonFamilies(aHit, cHit, gene).size() / 2) + "\t" + pvals.get(gene)));

		List<String> select = FDR.select(pvals, limits, 0.1);
		System.out.println("select.size() = " + select.size());

		return new Map[]{pvals, limits};
	}


	private static void printBiasOfAMutexGroup() throws IOException
	{
		String dir = "/home/ozgun/Documents/Grants/Mental/mutex/Groups/SFARI-10000/";
//		String dir = "/home/ozgun/Documents/Grants/Mental/mutex/Groups/Reactome-10000/";
//		String rID = "http://identifiers.org/reactome/R-HSA-400253"; // Circadian Clock
//		String rID = "http://identifiers.org/reactome/R-HSA-1257604"; // PI3K/AKT Signaling
		String rID = "1"; // PI3K/AKT Signaling
		String pFile = dir + escape(rID) + "-mutex.txt";

//		Set<String> genes = ReactomePathway.get().getGenes(rID);
		Set<String> genes = SFARI.get().getGenesWithMaxScore(1);

		Map<String, Set<String>> aHit = readHitMap(DenovoDB.DataFilterEnum.YUEN_TURNER_AUTISM);
		Map<String, Set<String>> cHit = readHitMap(DenovoDB.DataFilterEnum.YUEN_TURNER_CONTROL);

//		Map<String, Set<String>> aHit = readHitMap(DenovoDB.YUEN_AUTISM);
//		Map<String, Set<String>> cHit = readHitMap(DenovoDB.TURNER_AUTISM);

		int aTotal = getMutationCountInHitMap(aHit);
		int cTotal = getMutationCountInHitMap(cHit);

		double bias = cTotal / (double) (aTotal + cTotal);
		System.out.println("bias = " + bias);

		Map<String, Double> pMap = Files.lines(Paths.get(pFile)).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));

		List<String> ranked = new ArrayList<>(genes);
		Collections.shuffle(ranked);
		ranked = ranked.stream().sorted(Comparator.comparing(g -> pMap.getOrDefault(g, 1D))).collect(Collectors.toList());

		getMaxBiasOfRankedGeneList(ranked, aHit, cHit, bias);
	}

	private static String escape(String s)
	{
		return s.replaceAll(":", "_").replaceAll("/", "_");
	}

	private static void printBiasOfMutexGroups() throws IOException
	{
		String dir = "/home/ozgun/Documents/Temp/mssng-autism/Reactome-mssng-autism/";
//		String dir = "/home/ozgun/Documents/Temp/mssng-ssc-autism/SFARI-mssng-ssc-autism/";

		List<String[]> pathways = Files.lines(Paths.get(dir + "results.txt")).skip(1)
			.map(l -> l.split("\t")).collect(Collectors.toList());

//		String pathwayID = "R-HSA-400253"; // Circadian Clock
		String pathwayID = "R-HSA-1257604"; // PI3K/AKT Signaling
//		String pathwayID = "SFARI-1-to-1"; // SFARI rank 1
//		String pathwayID = "http://identifiers.org/reactome/R-HSA-983231"; // Platelet
//		String pathwayID = "http://identifiers.org/reactome/R-HSA-5620912"; // Anchoring

		Set<String[]> useSet = pathways.stream().filter(t -> t[0].equals(pathwayID)).collect(Collectors.toSet());

//		Set<String[]> useSet = pathways.stream()
//			.sorted((t1, t2) -> new Double(Double.valueOf(t2[3]) + Double.valueOf(t2[4])).compareTo(
//				Double.valueOf(t1[3]) + Double.valueOf(t1[4]))).limit(50).collect(Collectors.toSet());

		Map<String, Set<String>> xHit = readHitMap(DenovoDB.DataFilterEnum.YUEN_AUTISM);
//		Map<String, Set<String>> aHit = readHitMap(DenovoDB.DataFilterEnum.AN_AUTISM);
//		Map<String, Set<String>> cHit = readHitMap(DenovoDB.DataFilterEnum.AN_CONTROL);
		Map<String, Set<String>> aHit = readHitMap(DenovoDB.DataFilterEnum.TURNER_AUTISM);
		Map<String, Set<String>> cHit = readHitMap(DenovoDB.DataFilterEnum.AN_DIF_AUTISM);

//		Map<String, Set<String>> aHit = readHitMap(DenovoDB.YUEN_AUTISM);
//		Map<String, Set<String>> cHit = readHitMap(DenovoDB.TURNER_AUTISM);

		int aTotal = getMutationCountInHitMap(aHit);
		int cTotal = getMutationCountInHitMap(cHit);

		System.out.println("aTotal = " + aTotal);
		System.out.println("cTotal = " + cTotal);

		double bias = cTotal / (double) (aTotal + cTotal);
		System.out.println("bias = " + bias);

		Map<String, Double> pvals = new HashMap<>();

		for (String[] pathway : pathways)
		{
			if (!useSet.contains(pathway)) continue;

			String file = dir + pathway[0].replaceAll(":", "_").replaceAll("/", "_") + "-mutex.txt";

			List<String> geneList = Files.lines(Paths.get(file)).map(l -> l.split("\t")[0]).collect(Collectors.toList());

			if (pathway[0].startsWith("R-HSA"))
			{
				Set<String> genes = ReactomePathway.get().getGenes("http://identifiers.org/reactome/" + pathway[0]);
				genes.stream().filter(g -> !geneList.contains(g)).forEach(geneList::add);
			}

//			geneList.sort((g1, g2) -> Integer.compare(xHit.getOrDefault(g2, Collections.emptySet()).size(), xHit.getOrDefault(g1, Collections.emptySet()).size()));

			double p = getMaxBiasOfRankedGeneList(geneList, aHit, cHit, bias);
//			double p = getFinalBiasOfRankedGeneList(geneList, aHit, cHit, bias);
			pvals.put(pathway[0], p);

			System.out.println(pathway[5] + "\t" + p + "\t" + pathway[1]);
		}

		List<String> select = FDR.select(pvals, null, 0.1);
		System.out.println("select.size() = " + select.size());
//		Double thr = select.stream().map(pvals::get).max(Double::compareTo).get();
//		System.out.println("thr = " + thr);
	}

	private static double getMaxBiasOfRankedGeneList(List<String> geneList, Map<String, Set<String>> aHit, Map<String, Set<String>> cHit, double bias)
	{
		List<Double> pvals = new ArrayList<>();

		for (int i = 1; i <= geneList.size(); i++)
		{
			int aCnt = (int) geneList.stream().limit(i).map(g -> aHit.getOrDefault(g, Collections.emptySet()))
				.flatMap(Collection::stream).distinct().count();
			int cCnt = (int) geneList.stream().limit(i).map(g -> cHit.getOrDefault(g, Collections.emptySet()))
				.flatMap(Collection::stream).distinct().count();

			BinomialDistribution dist = new BinomialDistribution(aCnt + cCnt, bias);
			pvals.add(dist.cumulativeProbability(cCnt));
		}

		for (int i = 0; i < pvals.size(); i++) System.out.println((i + 1) + "\t" + pvals.get(i));

		return pvals.stream().min(Double::compareTo).get();
	}

	private static double getFinalBiasOfRankedGeneList(List<String> geneList, Map<String, Set<String>> aHit, Map<String, Set<String>> cHit, double bias)
	{
		int aCnt = (int) geneList.stream().map(g -> aHit.getOrDefault(g, Collections.emptySet()))
			.flatMap(Collection::stream).distinct().count();
		int cCnt = (int) geneList.stream().map(g -> cHit.getOrDefault(g, Collections.emptySet()))
			.flatMap(Collection::stream).distinct().count();

		BinomialDistribution dist = new BinomialDistribution(aCnt + cCnt, bias);
		return dist.cumulativeProbability(cCnt);
	}

	private static void printBiasOfReactomePathways() throws IOException
	{
		Map<String, Set<String>> aHit = readHitMap(DenovoDB.DataFilterEnum.YUEN_TURNER_AUTISM);
		Map<String, Set<String>> cHit = readHitMap(DenovoDB.DataFilterEnum.YUEN_TURNER_CONTROL);

//		keepPairedSamplesRemoveOverlaps(aHit, cHit);

		int aTotal = getMutationCountInHitMap(aHit);
		int cTotal = getMutationCountInHitMap(cHit);

		double bias = cTotal / (double) (aTotal + cTotal);
		System.out.println("bias = " + bias);

		Map<String, Set<String>> pathways = ReactomePathway.get().getAllPathways();

		Map<String, Double> pvals = new HashMap<>();
		Map<String, Double> limits = new HashMap<>();

		pathways.forEach((id, genes) ->
		{
			int aCnt = (int) genes.stream().map(g -> aHit.getOrDefault(g, Collections.emptySet()))
				.flatMap(Collection::stream).distinct().count();
			int cCnt = (int) genes.stream().map(g -> cHit.getOrDefault(g, Collections.emptySet()))
				.flatMap(Collection::stream).distinct().count();

			BinomialDistribution dist = new BinomialDistribution(aCnt + cCnt, bias);
			pvals.put(id, dist.cumulativeProbability(cCnt));
			limits.put(id, dist.cumulativeProbability(0));
		});

		List<String> select = FDR.select(pvals, limits, 0.1);
		select.forEach(id -> System.out.println( + pvals.get(id) + "\t" + id + "\t" + ReactomePathway.get().getName(id)));
	}

	private static void printDatasetBiasOfReactomePathways() throws IOException
	{
		DenovoDB.DataFilter filter = e -> DenovoDB.SELECT_MUT.select(e) && DenovoDB.SELECT_AUTISM_STUDY.select(e) &&
			e.studyName.equals("Yuen2017") && e.primaryPhenotype.equals("autism");

		Map<String, Set<String>> aHit = readHitMap(filter);

		filter = e -> DenovoDB.SELECT_MUT.select(e) && DenovoDB.SELECT_AUTISM_STUDY.select(e) &&
			e.studyName.equals("Turner_2017") && e.primaryPhenotype.equals("autism");

		Map<String, Set<String>> cHit = readHitMap(filter);

		int aTotal = getMutationCountInHitMap(aHit);
		int cTotal = getMutationCountInHitMap(cHit);

		double bias = cTotal / (double) (aTotal + cTotal);
		System.out.println("bias = " + bias);

		Map<String, Set<String>> pathways = ReactomePathway.get().getAllPathways();

		Map<String, Double> pvals = new HashMap<>();
		Map<String, Double> limits = new HashMap<>();

		pathways.forEach((id, genes) ->
		{
			int aCnt = (int) genes.stream().map(g -> aHit.getOrDefault(g, Collections.emptySet()))
				.flatMap(Collection::stream).distinct().count();
			int cCnt = (int) genes.stream().map(g -> cHit.getOrDefault(g, Collections.emptySet()))
				.flatMap(Collection::stream).distinct().count();

			BinomialDistribution dist = new BinomialDistribution(aCnt + cCnt, bias);
			pvals.put(id, dist.cumulativeProbability(cCnt));
			limits.put(id, dist.cumulativeProbability(0));
		});

		List<String> select = FDR.select(pvals, limits, 0.1);
		select.forEach(id -> System.out.println( + pvals.get(id) + "\t" + id + "\t" + ReactomePathway.get().getName(id)));
	}

	private static void compareUniformityVErsusFreqROCs() throws IOException
	{
		List<String> uniList = Files.lines(Paths.get(
			"/home/ozgun/Documents/Grants/Mental/mutation-heterogenetity-ranked.txt"))
			.map(l -> l.split("\t")[0]).collect(Collectors.toList());

		List<String> freqList = Files.lines(Paths.get(
			"/home/ozgun/Documents/Grants/Mental/number-of-mutations-ranked.txt"))
			.map(l -> l.split("\t")[0]).collect(Collectors.toList());

		Set<String> sfariGenes = new HashSet<>(SFARI.get().getGenesWithMaxScore(3));

		System.out.println("autism auc = " + ROC.getAUC(uniList, sfariGenes));
		System.out.println("control auc = " + ROC.getAUC(freqList, sfariGenes));

		PrintStream out = new PrintStream(new FileOutputStream("/home/ozgun/Documents/Grants/Mental/ROC2.txt"));
		ROC.printPlotForGoogleSheets(sfariGenes, out, uniList, freqList);
	}

	private static List<String> readMutexResultsAsRankedList(String file) throws IOException
	{
		List<String> ranked = new ArrayList<>();
		Set<String> seen = new HashSet<>();
		Files.lines(Paths.get(file)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			for (int i = 1; i < t.length; i++)
			{
				if (!seen.contains(t[i]))
				{
					seen.add(t[i]);
					ranked.add(t[i]);
				}
			}
		});
		return ranked;
	}

	private static void testMutexnessOfSFARI() throws IOException
	{
		List<String> ranks = Arrays.asList("1", "2", "3", "4", "5", "6", "all");
//		List<String> ranks = Arrays.asList("S", "1", "S1", "2", "S2", "3", "S3", "4", "S4", "5", "S5", "6", "S6", "all");

		String outDir = "/home/ozgun/Documents/Grants/Mental/mutex/Groups/SFARI-10000-Yuen-autism-Turner-control";
		FileUtil.mkdirs(outDir);

		Map<String, boolean[]> matrix = generateAlterationMatrixFromDenovoDB(null, DenovoDB.DataFilterEnum.YUEN_AUTISM_TURNER_CONTROL);

		Map<String, Integer> currentCovMap = new HashMap<>();
		Map<String, Integer> currentOvMap = new HashMap<>();
		Map<String, Set<String>> genesMap = new HashMap<>();

		ranks.forEach(rank ->
		{
			Set<String> genes = rank.equals("all") ? new HashSet<>(SFARI.get().getAllGenes()) :
				rank.equals("S") ? SFARI.get().getGenesSyndromic() :
					rank.startsWith("S") ? SFARI.get().getGenesWithMaxScore(Integer.valueOf(rank.substring(1))) :
						SFARI.get().getGenesWithMaxScore(Integer.valueOf(rank));

			if (rank.startsWith("S") && rank.length() > 1) genes.addAll(SFARI.get().getGenesSyndromic());

			genes.retainAll(matrix.keySet());
			genesMap.put(rank, genes);
		});

		// merge some gene clusters
		makeAdjustmentsToAlterationMatrixAndGeneSets(matrix, genesMap);

//		genesMap.putAll(AltMatrixUtil.getRandomControlsWithSimilarCoverages(matrix, genesMap, genesMap.get("all")));

		ranks = new ArrayList<>(genesMap.keySet());

		genesMap.forEach((rank, genes) ->
		{
			currentCovMap.put(rank, countCoverage(matrix, genes));
			currentOvMap.put(rank, countOverlap(matrix, genes));
		});

		Map<String, Double>[] pvals = AltMatrixUtil.getMutexCoocPvalsPreserveSampleWeights(matrix, genesMap, 10000, outDir);

		BufferedWriter writer = FileUtil.newBufferedWriter(outDir + "/results.txt");
		writer.write("SFARI Rank\tGenes size\tCoverage\tOverlap\tMutex p-value\tCooc p-value");
		ranks.stream().sorted(Comparator.comparing(n -> pvals[0].get(n))).forEach(rank -> FileUtil.lnwrite(
			rank + "\t" + genesMap.get(rank).size() + "\t" + currentCovMap.get(rank) + "\t" +
				currentOvMap.get(rank) + "\t" + pvals[0].get(rank) + "\t" + pvals[1].get(rank), writer));
		writer.close();
	}

	private static void testMutexnessOfReactomePathways() throws IOException
	{
		String outDir = "/home/ozgun/Documents/Grants/Mental/mutex/Groups/Reactome-10000-Yuen-autism-Turner-control";
		FileUtil.mkdirs(outDir);

		Map<String, boolean[]> matrix = generateAlterationMatrixFromDenovoDB(null, DenovoDB.DataFilterEnum.YUEN_AUTISM_TURNER_CONTROL);

		Map<String, Set<String>> genesMap = ReactomePathway.get().getCroppedPathways(matrix.keySet());
		System.out.println("genesMap.size() = " + genesMap.size());

		Map<String, Integer> currentCovMap = new HashMap<>();
		Map<String, Integer> currentOvMap = new HashMap<>();

		// merge some gene clusters
		makeAdjustmentsToAlterationMatrixAndGeneSets(matrix, genesMap);

		genesMap.forEach((id, genes) ->
		{
			currentCovMap.put(id, countCoverage(matrix, genes));
			currentOvMap.put(id, countOverlap(matrix, genes));
		});

		Map<String, Double>[] pvals = AltMatrixUtil.getMutexCoocPvalsPreserveSampleWeights(matrix, genesMap, 10000, outDir);

		BufferedWriter writer = FileUtil.newBufferedWriter(outDir + "/results.txt");
		writer.write("ID\tName\tGenes size\tCoverage\tOverlap\tMutex p-value\tCooc p-value");
		genesMap.keySet().stream().sorted(Comparator.comparing(n -> pvals[0].get(n))).forEach(id -> FileUtil.lnwrite(
			id + "\t" + ReactomePathway.get().getName(id) + "\t" + genesMap.get(id).size() + "\t" + currentCovMap.get(id) + "\t" +
				currentOvMap.get(id) + "\t" + pvals[0].get(id) + "\t" + pvals[1].get(id), writer));
		writer.close();
	}

	private static void makeAdjustmentsToAlterationMatrixAndGeneSets(Map<String, boolean[]> matrix,
		Map<String, Set<String>> genesMap)
	{
		Set<String> pcdhaGenes = FileUtil.getTermsInTabDelimitedColumn(
			"/home/ozgun/Documents/Grants/Mental/mutex/PCDH-cluster-genes.txt", 0, 0);

		boolean[] b = null;
		for (String gene : pcdhaGenes)
		{
			if (matrix.containsKey(gene))
			{
				if (b == null) b = matrix.get(gene);
				else ArrayUtil.ORWith(b, matrix.get(gene));
			}
		}
		if (b != null)
		{
			pcdhaGenes.forEach(matrix::remove);
			String pcdhaClusterGene = "PCDHA-cluster";
			matrix.put(pcdhaClusterGene, b);

			for (Set<String> set : genesMap.values())
			{
				if (set.removeAll(pcdhaGenes)) set.add(pcdhaClusterGene);
			}
		}

//		List<String> remGenes = Arrays.asList("PARD3B", "OTUD7A");
//		remGenes.forEach(matrix::remove);
//		genesMap.values().forEach(set -> set.removeAll(remGenes));
	}

	private static Map<String, boolean[]> generateAlterationMatrixFromDenovoDB(Set<String> genes, DenovoDB.DataFilter filter)
	{
		try
		{
			Map<String, Set<String>> hitMap = new HashMap<>();

			DenovoDB.get().getDataStream(filter).forEach(e ->
			{
				if (!hitMap.containsKey(e.gene)) hitMap.put(e.gene, new HashSet<>());
				hitMap.get(e.gene).add(e.sampleID);
			});

			BufferedWriter writer = Files.newBufferedWriter(Paths.get(
				"/home/ozgun/Documents/Grants/Mental/mutex/data-matrix.txt"));

			List<String> samples = DenovoDB.get().getDataStream(filter).map(e -> e.sampleID).distinct().sorted()
				.collect(Collectors.toList());

			System.out.println("samples.size() = " + samples.size());

			Map<String, boolean[]> matrix = new HashMap<>();

			for (String sample : samples)
			{
				writer.write("\t" + sample);
			}

			hitMap.keySet().stream().sorted((o1, o2) -> Integer.compare(hitMap.get(o2).size(), hitMap.get(o1).size()))
				.filter(g -> genes == null || genes.contains(g)).forEach(g ->
			{
				boolean[] b = new boolean[samples.size()];

				FileUtil.lnwrite(g, writer);

				for (int i = 0; i < samples.size(); i++)
				{
					b[i] = hitMap.get(g).contains(samples.get(i));
					FileUtil.tab_write(b[i] ? "1" : "0", writer);
				}
				matrix.put(g, b);
			});

			writer.close();
			return matrix;
		}
		catch (IOException e)
		{
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

	private static double getOverallMutexnessPValue(List<boolean[]> m)
	{
		int current = countCoverage(m);

		Histogram h = new Histogram(1);
		h.setBorderAtZero(true);
		h.setUseLowerBorderForPrinting(true);
		int trials = 10000;
		int meet = 0;
		for (int i = 0; i < trials; i++)
		{
			for (boolean[] b : m)
			{
				shuffleArray(b);
			}

			int c = countCoverage(m);
			h.count(c);
			if (c >= current) meet++;
		}
//		h.printDensity();

		return meet / (double) trials;
	}

	// Implementing Fisher–Yates shuffle
	static void shuffleArray(boolean[] b)
	{
		Random rnd = new Random();
		for (int i = b.length - 1; i > 0; i--)
		{
			int index = rnd.nextInt(i + 1);
			// Simple swap
			boolean a = b[index];
			b[index] = b[i];
			b[i] = a;
		}
	}

	static void randomizeMatrixPreserveSampleAlterations(List<boolean[]> m)
	{
		List<Integer> indices = new ArrayList<>();
		for (int i = 0; i < m.size(); i++)
		{
			indices.add(i);
		}

		Random r = new Random();

		for (boolean[] b1 : m)
		{
			for (int j = b1.length - 1; j > 0; j--)
			{
				int x = r.nextInt(j);

				if (b1[j] != b1[x])
				{
					Collections.shuffle(indices);

					boolean[] b2 = null;
					for (Integer indice : indices)
					{
						boolean[] b = m.get(indice);

						if (b1[j] == b[x] && b1[x] == b[j])
						{
							b2 = b;
							break;
						}
					}

					if (b2 != null)
					{
						b1[j] = !b1[j];
						b1[x] = !b1[x];
						b2[j] = !b2[j];
						b2[x] = !b2[x];
					}
				}
			}
		}
	}



	private static int countCoverage(List<boolean[]> matrix)
	{
		int n = matrix.get(0).length;
		int cov = 0;
		for (int i = 0; i < n; i++)
		{
			for (boolean[] b : matrix)
			{
				if (b[i])
				{
					cov++;
					break;
				}
			}
		}
		return cov;
	}

	private static int countCoverage(Map<String, boolean[]> matrix, Set<String> genes)
	{
		int n = matrix.get(genes.iterator().next()).length;
		int cov = 0;
		for (int i = 0; i < n; i++)
		{
			for (String gene : genes)
			{
				if (matrix.get(gene)[i])
				{
					cov++;
					break;
				}
			}
		}
		return cov;

	}

	private static int countOverlap(List<boolean[]> matrix)
	{
		int n = matrix.get(0).length;
		int ov = 0;

		for (int i = 0; i < n; i++)
		{
			boolean covered = false;
			for (boolean[] b : matrix)
			{
				if (b[i])
				{
					if (!covered) covered = true;
					else ov++;
				}
			}
		}
		return ov;
	}

	private static int countOverlap(Map<String, boolean[]> matrix, Set<String> genes)
	{
		int n = matrix.get(genes.iterator().next()).length;
		int ov = 0;

		for (int i = 0; i < n; i++)
		{
			boolean covered = false;
			for (String gene : genes)
			{
				if (matrix.get(gene)[i])
				{
					if (!covered) covered = true;
					else ov++;
				}
			}
		}
		return ov;
	}

	private static void printSignificantPairwiseCooc(Map<String, boolean[]> matrix)
	{
		System.out.println("Significant coocs:");
		Map<String, Double> pvals = new HashMap<>();

		for (String gene1 : matrix.keySet())
		{
			for (String gene2 : matrix.keySet())
			{
				if (gene1.compareTo(gene2) < 0)
				{
					double p = Overlap.calcCoocPval(matrix.get(gene1), matrix.get(gene2));
					pvals.put(gene1 + "   " + gene2, p);
				}
			}
		}

		List<String> select = FDR.select(pvals, null, 0.1);
		select.forEach(k -> System.out.println(k + "\t" + pvals.get(k)));
	}

	private static void printOverlapsOfMutations()
	{
		Map<String, Set<String>> hitMap = readHitMapForAutismFromDenovoDB();
//		Map<String, Set<String>> hitMap = readHitMapForControlsInAutismStudiesFromDenovoDB();

		Set<String> genes = SFARI.get().getGenesWithMaxScore(2);
//		Set<String> genes = hitMap.keySet();

		for (String gene1 : genes)
		{
			if (hitMap.containsKey(gene1))
			{
				Set<String> hits1 = hitMap.get(gene1);
				for (String gene2 : genes)
				{
					if (gene1.compareTo(gene2) < 0 && hitMap.containsKey(gene2))
					{
						Set<String> hits2 = hitMap.get(gene2);
						int o = CollectionUtil.countOverlap(hits1, hits2);

						if (o > 0)
						{
							System.out.println(gene1 + "\t" + gene2 + "\t" + o + "\t" + CollectionUtil.getIntersection(hits1, hits2));
						}
					}
				}
			}
		}
	}

	private static void selectFromMutexResults() throws IOException
	{
		for (int i = 0; i < 500; i+=10)
		{
			int hitThr = i;
			Map<String, Double> pvals = Files.lines(Paths.get(
				"/home/ozgun/Documents/Grants/Mental/mutex/Groups/Reactome-10000/results.txt"))
				.skip(1).map(l -> l.split("\t"))
				.filter(t -> ((Integer.valueOf(t[3]) + Integer.valueOf(t[4])) /*/ (double) Integer.valueOf(t[2])*/) > hitThr)
				.collect(Collectors.toMap(t -> t[1], t -> Double.valueOf(t[5]), Math::min));

			List<String> select = FDR.select(pvals, null, 0.1);
			System.out.println(i + "\t" + pvals.size() + "\t" + select.size());
		}
//		select.forEach(name -> System.out.println(pvals.get(name) + "\t" + name));
	}

	private static void printSelectPathways() throws IOException
	{

		String filename = "/home/ozgun/Documents/Grants/Mental/mutex/Groups/Reactome-10000/results.txt";
		int selectSize = 50;
		double fdrThr = 0.1;

		Map<String, Integer> idToHit = Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> Integer.valueOf(t[3]) + Integer.valueOf(t[4])));

		List<String> ranked = idToHit.keySet().stream().sorted(Comparator.comparing(idToHit::get).reversed()).collect(Collectors.toList());
		List<String> selected = ranked.subList(0, selectSize);

		Map<String, Double> pvals = Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t"))
			.filter(t -> selected.contains(t[0])).collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[5])));

		List<String> significant = FDR.select(pvals, null, fdrThr);
		System.out.println("significant.size() = " + significant.size());

		Files.lines(Paths.get(filename)).skip(1).forEach(l ->
		{
			String[] t = l.split("\t");
			if (selected.contains(t[0]))
			{
				System.out.println(l);
			}
		});
	}

	private static void printSFARIRanksOfGenes() throws IOException
	{
//		Map<String, Double> biasPvals = Files.lines(Paths.get("/home/ozgun/Documents/Grants/Mental/Turner-gene-bias.txt"))
//			.skip(1).map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[1], t -> Double.valueOf(t[4])));

		String groupFile = "/home/ozgun/Documents/Grants/Mental/mutex/Groups/SFARI-10000/" +
//			escape("http://identifiers.org/reactome/R-HSA-400253") // Circadian Clock
//			escape("http://identifiers.org/reactome/R-HSA-1257604") // PI3K/AKT Signaling
			escape("3") // rank 3
//		+ "-mutex.txt";
		+ "-cooc.txt";

		Set<String> genes = Files.lines(Paths.get(groupFile)).map(l -> l.split("\t")[0]).collect(Collectors.toSet());

		Map<String, boolean[]> matrix = generateAlterationMatrixFromDenovoDB(genes, DenovoDB.DataFilterEnum.YUEN_TURNER_AUTISM);

		Files.lines(Paths.get(groupFile)).map(l -> l.split("\t")).forEach(t ->
			System.out.println(SFARI.get().getClassification(t[0]) + "\t" + t[0] + "\t" + ArrayUtil.countValue(matrix.get(t[0]), true) + "\t" + t[1]));
//				+ "\t" + biasPvals.get(t[0])));
	}

	private static void findSignificantMembersOfGroups() throws IOException
	{
		String dir = "/home/ozgun/Documents/Grants/Mental/mutex/Groups/SFARI-10000/";

		for (File file : new File(dir).listFiles())
		{
			if (file.getName().endsWith("cooc.txt"))
			{
				Map<String, Double> pvals = Files.lines(Paths.get(file.getPath())).map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));
				List<String> select = FDR.select(pvals, null, 0.2);

				if (!select.isEmpty())
				{
					System.out.println("File = " + file.getName());
					select.forEach(g -> System.out.println(g + "\t" + pvals.get(g)));
				}
			}
		}
	}


	private static void compareTurnerAndAn() throws IOException
	{
		Set<String> samples = DenovoDB.get().getDataStream(e -> e.studyName.equals("Turner_2017"))
			.map(e -> e.sampleID).collect(Collectors.toSet());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get("/home/ozgun/Documents/Temp/turner-an-compare.txt"));
		writer.write("Sample\tCnt T only\tCnt Common\tCnt A only\tPresent in Turner but absent in An\tPresent in An but absent in Turner");

		for (String sample : samples)
		{
			Set<String> turnerGenes = DenovoDB.get().getDataStream(e -> e.sampleID.equals(sample) && e.studyName.equals("Turner_2017"))
				.map(e -> e.gene).filter(g -> !g.isEmpty() && !g.equals("NA")).collect(Collectors.toSet());
			Set<String> anGenes = DenovoDB.get().getDataStream(e -> e.sampleID.equals(sample) && e.studyName.equals("An2018"))
				.map(e -> e.gene).filter(g -> !g.isEmpty() && !g.equals("NA")).collect(Collectors.toSet());

			if (!turnerGenes.equals(anGenes))
			{
				List<String> diffTA = CollectionUtil.diff(turnerGenes, anGenes).stream().sorted().collect(Collectors.toList());
				List<String> diffAT = CollectionUtil.diff(anGenes, turnerGenes).stream().sorted().collect(Collectors.toList());
				writer.write("\n" + sample + "\t" + diffTA.size() + "\t" + CollectionUtil.countOverlap(anGenes, turnerGenes) +
					"\t" + diffAT.size() + "\t" + diffTA + "\t" + diffAT);
			}
		}

		writer.close();
	}
}
