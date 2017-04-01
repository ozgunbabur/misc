package org.panda.misc.pancan;

import org.panda.misc.MutexReader;
import org.panda.resource.CancerGeneBushman;
import org.panda.resource.CancerGeneCensus;
import org.panda.resource.OncoKB;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;
import org.panda.utility.statistics.FDR;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by babur on 2/23/16.
 */
public class MutexResultAggregator
{
	public static void main(String[] args) throws IOException
	{
//		printGroupsOf("/home/babur/Documents/PanCan/PanCan-results/", "TG", 0.1);

		int maxDiv = 30;
		prepareGeneRanking("/home/babur/Documents/PanCan/tissue-unnormalized-results/PanCan-results/", "/home/babur/Documents/PanCan/tissue-unnormalized-results/pancan-2.txt", maxDiv); System.out.println("---");
//		prepareGeneRanking("/home/babur/Documents/PanCan/PanCan-shuffled-1-results/", "/home/babur/Documents/PanCan/pancan-shuffled.txt", maxDiv);

//		checkSignalToNoise();

//		doReachAssessmentHardScoreThreshold();
	}

	private static void printGroupsOf(String dir, String gene, double scoreThr) throws IOException
	{
		Set<String> allGenes = new HashSet<>();
		Set<MutexReader.Group> results = MutexReader.readMutexResultsRecursive(dir, new HashSet<>());
		for (MutexReader.Group group : results)
		{
			if (group.genes.contains(gene) && group.score <= scoreThr)
			{
				System.out.println(group.toString().replaceAll(dir, ""));
				allGenes.addAll(group.genes);
			}
		}
		System.out.println();
		allGenes.stream().peek(System.out::print).forEach(g -> System.out.print(" "));
	}

	private static Map<String, Set<String>> getGeneParticipantsMerged(Map<String, Map<String, Double>> pairScores,
		double scoreThr) throws IOException
	{
		Map<String, Set<String>> map = new HashMap<>();

		pairScores.keySet().forEach(g1 ->
			map.put(g1, pairScores.get(g1).keySet().stream()
				.filter(g2 -> pairScores.get(g1).get(g2) <= scoreThr).collect(Collectors.toSet())));

		return map;
	}

	private static void prepareGeneRanking(String dir, String outFile, int divThr) throws IOException
	{
		Set<MutexReader.Group> results = MutexReader.readMutexResultsRecursive(dir, new HashSet<>(), f -> !hasNameOverThan(f.getName(), divThr));
		Set<String> canGen = new HashSet<>(CancerGeneCensus.get().getAllSymbols());
		canGen.addAll(OncoKB.get().getAllSymbols());

		Map<String, Double> scoreMap = new HashMap<>();
		Map<String, String> locMap = new HashMap<>();

		for (MutexReader.Group group : results)
		{
			for (String gene : group.genes)
			{
				if (!scoreMap.containsKey(gene) || scoreMap.get(gene) > group.score)
				{
					scoreMap.put(gene, group.score);
					locMap.put(gene, group.fromDir.substring(dir.length()));
				}
			}
		}

//		Map<String, Double> secondBest = getSecondBestScoresOnNonOverlappingParts(results);
//		Map<String, Double> secondBest = getBestScores(results);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("Gene\tIn CGC\tScore\tChunk\tNumberOfNeighbors\tNames of neighbors");

//		Map<String, Double> pairScores = getSecondBestPairScoresOnNonOverlappingParts(results);
//		Map<String, Double> pairScores = getBestPairScores(results);

		scoreMap.keySet().stream().sorted((g1, g2) -> scoreMap.get(g1).compareTo(scoreMap.get(g2))).forEach(g ->
			FileUtil.lnwrite(g + "\t" + (canGen.contains(g) ? "X" : "") + "\t" +
				scoreMap.get(g) + "\t" + locMap.get(g), writer));// + "\t" +
//				getParticipants(g, pairScores, 0.01).size() + "\t" + getParticipants(g, pairScores, 0.01), writer));

		writer.close();

//		BufferedWriter writer2 = new BufferedWriter(new FileWriter(outFile.substring(0, outFile.lastIndexOf(".")) +
//			"-pair-scores.txt"));
//		pairScores.keySet().forEach(k -> FileUtil.writeln(k + "\t" + pairScores.get(k), writer2));
//		writer2.close();
	}

	private static List<String> getParticipants(String gene, Map<String, Double> pairScores, double thr)
	{
		return  pairScores.keySet().stream()
			.filter(k -> k.startsWith(gene + "\t") || k.endsWith("\t" + gene))
			.filter(k -> pairScores.get(k) <= thr)
			.map(k -> k.split("\t")).flatMap(Arrays::stream).filter(g -> !g.equals(gene))
			.sorted().collect(Collectors.toList());
	}

	private static String getParticipants(String gene, Map<String, Set<String>> map)
	{
		if (map.containsKey(gene))
		{
			return map.get(gene).toString().replaceAll(",", "").replaceAll("\\[", "").replaceAll("]", "");
		}
		return "";
	}

	private static Map<String, Double> getBestScores(Set<MutexReader.Group> results)
	{
		Map<String, Double> scores = new HashMap<>();

		for (MutexReader.Group group : results)
		{
			for (String gene : group.genes)
			{
				if (!scores.containsKey(gene) || scores.get(gene) > group.score) scores.put(gene, group.score);
			}
		}

		return scores;
	}

	private static List<Map<String, Double>> getBestScores(List<Set<MutexReader.Group>> resultsList)
	{
		return resultsList.stream().map(MutexResultAggregator::getBestScores).collect(Collectors.toList());
	}

	private static Map<String, Double> getSecondBestScoresOnNonOverlappingParts(Set<MutexReader.Group> results)
	{
		Map<String, Set<MutexReader.Group>> geneToGroups = new HashMap<>();
		for (MutexReader.Group group : results)
		{
			for (String gene : group.genes)
			{
				if (!geneToGroups.containsKey(gene)) geneToGroups.put(gene, new HashSet<>());
				geneToGroups.get(gene).add(group);
			}
		}

		return getNonOverlappingSecondBestScores(geneToGroups);
	}

	private static Map<String, Double> getBestPairScores(Set<MutexReader.Group> results)
	{
		Map<String, Set<MutexReader.Group>> pairToGroups = getPairToGroupsMap(results);
		return getBestScores(pairToGroups);
	}

	private static Map<String, Double> getSecondBestPairScoresOnNonOverlappingParts(Set<MutexReader.Group> results)
	{
		Map<String, Set<MutexReader.Group>> pairToGroups = getPairToGroupsMap(results);
		return getNonOverlappingSecondBestScores(pairToGroups);
	}

	private static Map<String, Set<MutexReader.Group>> getPairToGroupsMap(Set<MutexReader.Group> results)
	{
		Map<String, Set<MutexReader.Group>> pairToGroups = new HashMap<>();
		for (MutexReader.Group group : results)
		{
			for (String g1 : group.genes)
			{
				for (String g2 : group.genes)
				{
					if (g1.compareTo(g2) < 0)
					{
						String key = g1 + "\t" + g2;
						if (!pairToGroups.containsKey(key)) pairToGroups.put(key, new HashSet<>());
						pairToGroups.get(key).add(group);
					}
				}
			}
		}
		return pairToGroups;
	}

	private static Map<String, Double> getNonOverlappingSecondBestScores(
		Map<String, Set<MutexReader.Group>> keyToGroups)
	{
		Progress p = new Progress(keyToGroups.size(), "Calculating second best scores");
		Map<String, Double> secondBest = new HashMap<>();
		for (String key : keyToGroups.keySet())
		{
			p.tick();
			Set<MutexReader.Group> groups = keyToGroups.get(key);
			for (MutexReader.Group g1 : groups)
			{
				if (secondBest.containsKey(key) && secondBest.get(key) <= g1.score) continue;

				for (MutexReader.Group g2 : groups)
				{
					if (secondBest.containsKey(key) && secondBest.get(key) <= g2.score) continue;
					if (g1.hashCode() <= g2.hashCode() && overlapFree(g1, g2))
					{
						double val = Math.max(g1.score, g2.score);
						if (!secondBest.containsKey(key) || secondBest.get(key) > val)
						{
							secondBest.put(key, val);
						}
					}
				}
			}
		}
		return secondBest;
	}

	private static Map<String, Double> getBestScores(Map<String, Set<MutexReader.Group>> keyToGroups)
	{
		Map<String, Double> best = new HashMap<>();
		for (String key : keyToGroups.keySet())
		{
			Set<MutexReader.Group> groups = keyToGroups.get(key);
			for (MutexReader.Group group : groups)
			{
				if (!best.containsKey(key) || best.get(key) > group.score) best.put(key, group.score);
			}
		}
		return best;
	}

	private static boolean overlapFree(MutexReader.Group g1, MutexReader.Group g2)
	{
		String[] split1 = g1.fromDir.split("/");
		String[] split2 = g2.fromDir.split("/");
		return overlapFree(Integer.parseInt(split1[6]), Integer.parseInt(split1[7]),
			Integer.parseInt(split2[6]), Integer.parseInt(split2[7]));
	}

	private static boolean overlapFree(int totalParts1, int part1, int totalParts2, int part2)
	{
		double[] borders1 = getBorders(totalParts1, part1);
		double[] borders2 = getBorders(totalParts2, part2);
		return borders1[1] < borders2[0] || borders1[0] > borders2[1];
	}

	private static double[] getBorders(int totalPats, int part)
	{
		double piece = 1D / totalPats;
		return new double[]{piece * (part - 1), piece * part};
	}

	private static void checkSignalToNoise()
	{
		double[] fdr = new double[20];
		for (int i = 0; i < 20; i++)
		{
			fdr[i] = (i + 1) * 0.01;
		}

		int[][] sizes = new int[30][fdr.length];

		String resultDir = "/home/babur/Documents/PanCan/PanCan-results";
		String shuffledDir = "/home/babur/Documents/PanCan/PanCan-shuffled-?-results";

		for (int i = 0; i < 30; i++)
		{
			int divisionThr = i + 1;
			Set<MutexReader.Group> test = MutexReader.readMutexResultsRecursive(resultDir, new HashSet<>(), f -> !hasNameOverThan(f.getName(), divisionThr));
			Map<String, Double> testScores = getBestScores(test);

			List<Set<MutexReader.Group>> ctrl = readShuffledResults(shuffledDir, f -> !hasNameOverThan(f.getName(), divisionThr));
			List<Map<String, Double>> shuffledScores = getBestScores(ctrl);

			List<Double> noiseList = new ArrayList<>();
			for (Map<String, Double> scores : shuffledScores)
			{
				noiseList.addAll(scores.values());
			}

			for (int j = 0; j < fdr.length; j++)
			{
				List<String> select = FDR.select(testScores, fdr[j], noiseList, shuffledScores.size());
				sizes[i][j] = select.size();
			}
		}

		// normalize the matrix
		double[][] norm = new double[sizes.length][sizes[0].length];
		int[] max = new int[fdr.length];
		for (int j = 0; j < fdr.length; j++)
		{
			for (int i = 0; i < sizes.length; i++)
			{
				if (max[j] < sizes[i][j]) max[j] = sizes[i][j];
			}
		}
		for (int i = 0; i < sizes.length; i++)
		{
			for (int j = 0; j < sizes[i].length; j++)
			{
				norm[i][j] = sizes[i][j] / (double) max[j];
			}
		}

		System.out.print("Division");
		for (double aFdr : fdr)
		{
			System.out.print("\t" + aFdr);
		}

		for (int i = 0; i < sizes.length; i++)
		{
			System.out.print("\n" + (i + 1));
			for (int j = 0; j < fdr.length; j++)
			{
				System.out.print("\t" + norm[i][j]);
			}
		}
	}

	private static boolean hasNameOverThan(String name, int thr)
	{
		int x;
		try
		{
			x = Integer.parseInt(name);
		}
		catch (NumberFormatException e){return false;}

		return x > thr;
	}

	private static void doReachAssessmentHardScoreThreshold()
	{
		String testDir = "/home/babur/Documents/PanCan/PanCan-results";
		String ctrlDir = "/home/babur/Documents/PanCan/PanCan-shuffled-?-results";

		int[] noNetworkSize = new int[30];
		int[] withPCSize = new int[30];
		int[] allSize = new int[30];

		double[] noNetworkFDR = new double[30];
		double[] withPCFDR = new double[30];
		double[] allFDR = new double[30];

		double thr = 0.01;
		System.out.println("thr = " + thr);

		for (int i = 0; i < 30; i++)
		{
			int div = i + 1;

			Set<MutexReader.Group> noNetworkTestResults = MutexReader.readMutexResultsRecursive(testDir, new HashSet<>(), f -> !hasNameOverThan(f.getName(), div) && !f.getName().endsWith("PC2v8"));
			Map<String, Double> noNetworkTestScores = getBestScores(noNetworkTestResults);
			List<Set<MutexReader.Group>> noNetworkCtrlResults = readShuffledResults(ctrlDir, f -> !hasNameOverThan(f.getName(), div) && !f.getName().endsWith("PC2v8"));
			List<Map<String, Double>> noNetworkCtrlScores = getBestScores(noNetworkCtrlResults);

			noNetworkSize[i] = (int) noNetworkTestScores.keySet().stream().filter(g -> noNetworkTestScores.get(g) <= thr).count();
			double shufSize = getMeanShuffledSize(noNetworkCtrlScores, thr);
			noNetworkFDR[i] = shufSize / noNetworkSize[i];

			Set<MutexReader.Group> withPCTestResults = MutexReader.readMutexResultsRecursive(testDir, new HashSet<>(), f -> !hasNameOverThan(f.getName(), div) && !f.getName().startsWith("fries"));
			Map<String, Double> withPCTestScores = getBestScores(withPCTestResults);
			List<Set<MutexReader.Group>> withPCCtrlResults = readShuffledResults(ctrlDir, f -> !hasNameOverThan(f.getName(), div) && !f.getName().startsWith("fries"));
			List<Map<String, Double>> withPCCtrlScores = getBestScores(withPCCtrlResults);

			withPCSize[i] = (int) withPCTestScores.keySet().stream().filter(g -> withPCTestScores.get(g) <= thr).count();
			shufSize = getMeanShuffledSize(withPCCtrlScores, thr);
			withPCFDR[i] = shufSize / withPCSize[i];

			Set<MutexReader.Group> allTestResults = MutexReader.readMutexResultsRecursive(testDir, new HashSet<>(), f -> !hasNameOverThan(f.getName(), div));
			Map<String, Double> allTestScores = getBestScores(allTestResults);
			List<Set<MutexReader.Group>> allCtrlResults = readShuffledResults(ctrlDir, f -> !hasNameOverThan(f.getName(), div));
			List<Map<String, Double>> allCtrlScores = getBestScores(allCtrlResults);

			allSize[i] = (int) allTestScores.keySet().stream().filter(g -> allTestScores.get(g) <= thr).count();
			shufSize = getMeanShuffledSize(allCtrlScores, thr);
			allFDR[i] = shufSize / allSize[i];
		}

		System.out.print("Div\tNo-network size\tNo-network FDR\tWith-PC Size\tWith-PC FDR\tAll size\tAll FDR");
		for (int i = 0; i < 30; i++)
		{
			System.out.print("\n" + (i + 1));
			System.out.print("\t" + noNetworkSize[i] + "\t" + noNetworkFDR[i]);
			System.out.print("\t" + withPCSize[i] + "\t" + withPCFDR[i]);
			System.out.print("\t" + allSize[i] + "\t" + allFDR[i]);
		}
		System.out.println();
	}

	private static List<Set<MutexReader.Group>> readShuffledResults(String dirRoot, MutexReader.DirectoryFilter filter)
	{
		List<Set<MutexReader.Group>> list = new ArrayList<>();
		for (int i = 1; i < 100; i++)
		{
			String dir = dirRoot.replace("?", i + "");
			if (new File(dir).exists())
			{
				list.add(MutexReader.readMutexResultsRecursive(dir, new HashSet<>(), filter));
			}
		}
		return list;
	}

	private static double getMeanShuffledSize(List<Map<String, Double>> scoreList, double thr)
	{
		double total = 0;
		for (Map<String, Double> scores : scoreList)
		{
			total += scores.keySet().stream().filter(g -> scores.get(g) <= thr).count();
		}
		return total / scoreList.size();
	}


}
