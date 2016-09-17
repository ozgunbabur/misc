package org.panda.misc.pancan;

import org.panda.misc.MutexReader;
import org.panda.resource.CancerGeneBushman;
import org.panda.resource.CancerGeneCensus;
import org.panda.utility.*;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.Overlap;
import org.panda.utility.statistics.Summary;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class PanCanCooccurrenceTester
{
	public static final String COOC_FILENAME = "cooc-groups.txt";
	public static void main(String[] args) throws IOException
	{
//		double thr = 0.01;
//		writeCoocPvals("/home/babur/Documents/mutex/TCGA/PanCan-shuffled/20", thr, "/home/babur/Documents/PanCan/cooc-20-shuffled.txt");
//		compare();
//		writeROC();
//		printCosmicPresence();
//		prepareCoocSIF("/home/babur/Documents/PanCan/cooc-20.txt", 200, "/home/babur/Documents/PanCan/cooc-20-first200");
//		printDegreeRanking("/home/babur/Documents/PanCan/cooc-20.txt", 200);

//		Kronometre kron = new Kronometre();
//
//		String dir = "/home/babur/Documents/mutex/TCGA/PanCan";
//		writeCoocGroupsInTheirDirectoryRecursive(dir, dir, "/home/babur/Documents/PanCan/PanCan-results", 1000);
//
//		for (int i = 1; i <= 10; i++)
//		{
//			dir = "/home/babur/Documents/mutex/TCGA/PanCan-shuffled-" + i;
//			String destDir = "/home/babur/Documents/PanCan/PanCan-shuffled-" + i + "-results";
//
//			if (Files.exists(Paths.get(destDir)))
//			{
//				writeCoocGroupsInTheirDirectoryRecursive(dir, dir, destDir, 1000);
//			}
//		}
//		kron.stop();
//		kron.print();

//		resulSizeComparisons();

//		String base = "/home/babur/Documents/PanCan/";
//		Object[] o = readCoocResults(true, 30);
//		Set<MutexReader.Group> groups = (Set<MutexReader.Group>) o[0];
//		List<Double> ctrl = (List<Double>) o[1];
//		int randIter = (int) o[2];
//		writeRankedGenes(groups, ctrl, randIter, base + "cooc-genes-right.txt");
//		prepareCoocSIF(groups, ctrl, randIter, base + "cooc-right");
//		o = readCoocResults(false, 30);
//		groups = (Set<MutexReader.Group>) o[0];
//		ctrl = (List<Double>) o[1];
//		randIter = (int) o[2];
//		writeRankedGenes(groups, ctrl, randIter, base + "cooc-genes-left.txt");
//		prepareCoocSIF(groups, ctrl, randIter, base + "cooc-left");

//		printPerformanceForMaxDiv();

		printDiffFromMutexResults();
	}

	static void writeCoocGroupsInTheirDirectoryRecursive(String dir, String sourceBase, String destinationBase, int lineLimit) { try
	{
		if (Files.exists(Paths.get(dir + "/DataMatrix.txt")))
		{
			if (!Files.exists(Paths.get(dir.replace(sourceBase, destinationBase) + "/" + COOC_FILENAME)))
			{
				TopOfMap<String, Double> pvals = getPairCoocPVals(dir, lineLimit);
				BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir.replace(sourceBase, destinationBase) + "/" + COOC_FILENAME));
				writer.write("Score\tMembers");
				pvals.keySet().stream().sorted((k1, k2) -> pvals.get(k1).compareTo(pvals.get(k2)))
					.forEach(k -> FileUtil.lnwrite(pvals.get(k) + "\t" + k, writer));
				writer.close();
			}
		}
		else
		{
			Files.newDirectoryStream(Paths.get(dir)).forEach(p ->
			{
				if (Files.isDirectory(p))
				{
					writeCoocGroupsInTheirDirectoryRecursive(p.toString(), sourceBase, destinationBase, lineLimit);
				}
			});
		}
	} catch (IOException e){e.printStackTrace();}}


	static TopOfMap<String, Double> getPairCoocPVals(String leafDir, int sizeLimit) throws IOException
	{
		TopOfMap<String, Double> pvals = new TopOfMap<>(Double::compare, sizeLimit);

		System.out.print("Processing " + leafDir + " ... ");
		Map<String, boolean[]> altMap = readAlterations(leafDir + "/DataMatrix.txt");
		for (String g1 : altMap.keySet())
		{
			for (String g2 : altMap.keySet())
			{
				if (g1.compareTo(g2) < 0)
				{
					double pval = Overlap.calcCoocPval(altMap.get(g1), altMap.get(g2));
//					double pval = FishersExactTest.calcCoocPval(altMap.get(g1), altMap.get(g2));

					String key = g1 + "\t" + g2;
					pvals.put(key, pval);
				}
			}
		}
		System.out.println("ok");
		return pvals;
	}

	static Map<String, boolean[]> readAlterations(String file) throws IOException
	{
		Map<String, boolean[]> alts = new HashMap<>();
		Files.lines(Paths.get(file)).skip(1).map(l -> l.split("\t")).forEach(t -> alts.put(t[0], convertAlterations(t)));
		return alts;
	}

	static boolean[] convertAlterations(String[] token)
	{
		boolean[] b = new boolean[token.length - 1];
		for (int i = 1; i < token.length; i++)
		{
			b[i - 1] = !token[i].equals("0");
		}
		return b;
	}

	static void prepareCoocSIF(Set<MutexReader.Group> groups, List<Double> ctrl, int randIter, String sifWithoutExtension) throws IOException
	{
		Map<String, Double> geneScores = MutexReader.convertGroupsToGeneBestScores(groups);
		Map<String, Double> qVals = FDR.getQVals(geneScores, ctrl, randIter);

		Map<String, Double> edgeScores = new HashMap<>();
		groups.stream().forEach(g ->
		{
			String key = g.genes.get(0) + " in-same-group " + g.genes.get(1);
			if (!edgeScores.containsKey(key) || edgeScores.get(key) > g.getScore()) edgeScores.put(key, g.getScore());
		});

		double maxScore = geneScores.keySet().stream().filter(g -> qVals.get(g) <= 0.05).map(geneScores::get)
			.max(Double::compare).get();

		BufferedWriter writer1 = Files.newBufferedWriter(Paths.get(sifWithoutExtension + ".sif"));
		edgeScores.keySet().stream().filter(k -> edgeScores.get(k) <= maxScore).forEach(k ->
			FileUtil.writeln(k.replaceAll(" ", "\t"), writer1));
		writer1.close();

		Set<String> genesInGraph = geneScores.keySet().stream().filter(g -> qVals.get(g) <= 0.05).collect(Collectors.toSet());

		Map<String, Double> coverages = MutexToGraphAligner.readMutationCoverages(genesInGraph);

		ValToColor covVTC = new ValToColor(new double[]{0, 0.2}, new Color[]{Color.WHITE, Color.RED});

		double minScore = Summary.minOfDoubleCollection(geneScores.values());

		System.out.println("minScore = " + minScore);
		System.out.println("maxScore = " + maxScore);

		ValToColor signifVTC = new ValToColor(new double[]{-Math.log(maxScore) - 1, -Math.log(1E-4)}, new Color[]{Color.WHITE, Color.BLACK});
		NumberFormat fmt = new DecimalFormat("0.####");

		BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(sifWithoutExtension + ".format"));
		writer2.write("node\tall-nodes\tcolor\t255 255 255");
		writer2.write("\nnode\tall-nodes\tborderwidth\t2");
		writer2.write("\nedge\tall-edges\twidth\t2");
		for (String gene : genesInGraph)
		{
			writer2.write("\nnode\t" + gene + "\tcolor\t" + covVTC.getColorInString(coverages.get(gene)));
			writer2.write("\nnode\t" + gene + "\tbordercolor\t" + signifVTC.getColorInString(-Math.log(geneScores.get(gene))));
			writer2.write("\nnode\t" + gene + "\ttooltip\tp=" + fmt.format(geneScores.get(gene)) + ", c=" + fmt.format(coverages.get(gene)));
		}
		for (String edge : edgeScores.keySet())
		{
			writer2.write("\nedge\t" + edge + "\tcolor\t" + signifVTC.getColorInString(-Math.log(edgeScores.get(edge))));
		}
		writer2.close();
	}

	static void printDegreeRanking(String file, int lineLimit) throws IOException
	{
		Map<String, Set<String>> neighMap = new HashMap<>();
		Files.lines(Paths.get(file)).limit(lineLimit).map(l -> l.split("\t")).forEach(t ->
		{
			if (!neighMap.containsKey(t[0])) neighMap.put(t[0], new HashSet<>());
			neighMap.get(t[0]).add(t[1]);
			if (!neighMap.containsKey(t[1])) neighMap.put(t[1], new HashSet<>());
			neighMap.get(t[1]).add(t[0]);
		});

		neighMap.keySet().stream()
			.sorted((g1, g2) -> new Integer(neighMap.get(g2).size()).compareTo(neighMap.get(g1).size()))
			.forEach(g -> System.out.println(g + "\t" + neighMap.get(g).size() + "\t" + neighMap.get(g)));

		System.out.println("--- # of ZNFs = " + neighMap.keySet().stream().filter(g -> g.startsWith("ZNF")).count());
	}

	/**
	 * Object 1: Set<MutexReader.Group> actual results
	 * Object 2: List<Double> randomization values
	 * Object 3: Integer number of randomizations
	 */
	static Object[] readCoocResults(Boolean right, int maxDiv) throws IOException
	{
		String base = "/home/babur/Documents/PanCan/";
		String testDir = base + "PanCan-results";
		String ctrlDir = base + "PanCan-shuffled-?-results";

		MutexReader.DirectoryFilter filter = right == null ? f -> !violatesMaxDiv(f, maxDiv) :
			right ? f -> !violatesMaxDiv(f, maxDiv) && !isLeft(f) : f -> !violatesMaxDiv(f, maxDiv) && !isRight(f);

		return PanCanResultLoader.readGroupsWithFlattenedControl(false, testDir, ctrlDir, filter);
	}

	static boolean violatesMaxDiv(File f, int maxDiv)
	{
		int[] n = getBranchNums(f);
		return n != null && n[0] > maxDiv;
	}

	static boolean isRight(File f)
	{
		int[] n = getBranchNums(f);
		return n != null && n[0] == n[1];
	}

	static boolean isLeft(File f)
	{
		int[] n = getBranchNums(f);
		return n != null && n[0] != n[1];
	}

	static int[] getBranchNums(File f)
	{
		String[] t = f.getPath().split(File.separator);
		int n1;
		int n2;

		try
		{
			n2 = Integer.parseInt(t[t.length - 1]);
			n1 = Integer.parseInt(t[t.length - 2]);
		}
		catch (NumberFormatException e){ return null; }

		return new int[]{n1, n2};
	}

	private static void printPerformanceForMaxDiv() throws IOException
	{
		boolean roc = true;
		Boolean right = true;

		Map<String, List<Double>> testScoreMap = new LinkedHashMap<>();
		Map<String, List<Double>> ctrlScoreMap = new LinkedHashMap<>();
		double ss = 1;

		for (int i = 3; i < 30; i++)
		{
			int div = i + 1;

			Object[] o = readCoocResults(right, div);
			Set<MutexReader.Group> testResults = (Set<MutexReader.Group>) o[0];
			Map<String, Double> testScores = MutexReader.convertGroupsToGeneBestScores(testResults);
			List<Double> ctrlScores = (List<Double>) o[1];
			ss = (int) o[2];
			List<Double> testScoresList = new ArrayList<>(testScores.values());
			Collections.sort(testScoresList);
			Collections.sort(ctrlScores);
			testScoreMap.put(div + "", testScoresList);
			ctrlScoreMap.put(div + "", ctrlScores);
		}

		Map<String, Integer> cIndMap = new HashMap<>();
		for (String type : testScoreMap.keySet())
		{
			cIndMap.put(type, -1);
		}

		if (!roc) System.out.print("Size");

		for (String type : testScoreMap.keySet())
		{
			if (roc) System.out.print("\t" + type + "\t");
			else     System.out.print("\t" + type);
		}

		for (int j = 1; j <= 1000; j++)
		{
			System.out.println();
			if (!roc) System.out.print(j);

			for (String type : testScoreMap.keySet())
			{
				while (cIndMap.get(type) < ctrlScoreMap.get(type).size() - 1 &&
					ctrlScoreMap.get(type).get(cIndMap.get(type) + 1) < testScoreMap.get(type).get(j-1))
				{
					cIndMap.put(type, cIndMap.get(type) + 1);
				}

				if (roc)
				{
					double fp = (cIndMap.get(type) + 1) / ss;
					System.out.print(fp + "\t" + (j - fp) + "\t");
				}
				else
				{
					System.out.print("\t" + (((cIndMap.get(type) + 1) / ss) / j));
				}
			}
		}
	}



	static void writeRankedGenes(Set<MutexReader.Group> groups, List<Double> ctrl, int randIter, String filename) throws IOException
	{
		String common = "/home/babur/Documents/PanCan/PanCan-results/";
		Map<String, Set<String>> partners = new HashMap<>();
		Map<String, Double> geneScores = new HashMap<>();
		Map<String, String> geneLocations = new HashMap<>();
		groups.stream().peek(g -> g.shortenLoc(common)).forEach(g ->
		{
			for (String gene : g.genes)
			{
				if (!partners.containsKey(gene)) partners.put(gene, new HashSet<>());
				partners.get(gene).addAll(g.genes);

				if (!geneScores.containsKey(gene) || geneScores.get(gene) > g.getScore())
				{
					geneScores.put(gene, g.getScore());
					geneLocations.put(gene, g.fromDir);
				}
			}
		});

		Set<String> cancerGenes = new HashSet<>();
		cancerGenes.addAll(CancerGeneBushman.get().getAllSymbols());
		cancerGenes.addAll(CancerGeneCensus.get().getAllSymbols());

		Map<String, ArrayList<String>> partnerLists = partners.keySet().stream().peek(g -> partners.get(g).remove(g))
			.collect(Collectors.toMap(g -> g, g -> new ArrayList<>(partners.get(g))));

		partnerLists.keySet().forEach(g -> Collections.sort(partnerLists.get(g)));

		Map<String, Double> qVals = FDR.getQVals(geneScores, ctrl, randIter);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));
		writer.write("Gene\tIn CGC\tScore\tq-val\tBest Location\tDegree\tPartners");
		geneScores.keySet().stream().sorted((g1, g2) -> geneScores.get(g1).compareTo(geneScores.get(g2))).forEach(g ->
			FileUtil.lnwrite(g + "\t" + (cancerGenes.contains(g) ? "X" : "") + "\t" + geneScores.get(g) +
				"\t" + qVals.get(g) + "\t" + geneLocations.get(g) + "\t" + partners.get(g).size() + "\t" +
				partnerLists.get(g), writer));
		writer.close();
	}

	static void resulSizeComparisons() throws IOException
	{
		String base = "/home/babur/Documents/PanCan/";
		String testDir = base + "PanCan-results/";
		String ctrlDir = base + "PanCan-shuffled-1-results/";

		Map<String, Integer> cntMap = new HashMap<>();
		Map<String, Set<MutexReader.Group>> groupsMap = new HashMap<>();
		List<String> empty = new ArrayList<>();
		for (int i = 1; i <= 30; i++)
		{
			for (int j = 1; j <= i; j++)
			{
				String branch = i + "/" + j;
				Map<MutexReader.Group, Double> testScr = MutexReader.readCoocResults(testDir + branch).stream().collect(Collectors.toMap(g -> g, g -> g.score));
				List<Double> ctrlScr = Files.lines(Paths.get(ctrlDir + branch + "/cooc-groups.txt")).skip(1).filter(l -> !l.isEmpty()).map(l -> Double.parseDouble(l.split("\t")[0])).collect(Collectors.toList());
				Collections.sort(ctrlScr);

				List<MutexReader.Group> select = FDR.select(testScr, 0.001, ctrlScr, 1);
				if (!select.isEmpty())
				{
					cntMap.put(branch, select.size());
					groupsMap.put(branch, new HashSet<>(select));
				}
				else empty.add(branch);
			}
		}

		List<String> remove = new ArrayList<>();
		for (String b1 : groupsMap.keySet())
		{
			for (String b2 : groupsMap.keySet())
			{
				if (b1.equals(b2)) continue;

				Set<MutexReader.Group> g1 = groupsMap.get(b1);
				Set<MutexReader.Group> g2 = groupsMap.get(b2);
				if (g2.size() > g1.size() && g2.containsAll(g1)) remove.add(b1);
			}
		}

		cntMap.keySet().stream().sorted((b1, b2) -> cntMap.get(b2).compareTo(cntMap.get(b1))).filter(b -> !remove.contains(b)).forEach(b -> System.out.println(cntMap.get(b) + "\t" + b));

		System.out.println("empty = " + empty);
	}

	static void printDiffFromMutexResults() throws IOException
	{
		Set<String> mutex = Files.lines(Paths.get("/home/babur/Documents/PanCan/pancan.txt")).skip(1).limit(500)
			.map(l -> l.split("\t")[0]).collect(Collectors.toSet());

		Set<String> left = Files.lines(Paths.get("/home/babur/Documents/PanCan/cooc-genes-left.txt")).skip(1)
			.map(l -> l.split("\t")).filter(t -> Double.parseDouble(t[2]) <= 0.05).map(t -> t[0])
			.collect(Collectors.toSet());

		Set<String> right = Files.lines(Paths.get("/home/babur/Documents/PanCan/cooc-genes-right.txt")).skip(1)
			.map(l -> l.split("\t")).filter(t -> Double.parseDouble(t[2]) <= 0.05).map(t -> t[0])
			.collect(Collectors.toSet());

		CollectionUtil.printNameMapping("Mutex", "Cooc Left", "Cooc Right");
		CollectionUtil.printVennSets(mutex, left, right);
	}

	static void reportModifiedRands() throws FileNotFoundException
	{
		String base = "/home/babur/Documents/PanCan/PanCan-shuffled-1-results/";
		for (int i = 1; i <= 30; i++)
		{
			for (int j = 1; j <= i; j++)
			{
				String branch = i + "/" + j;
				Scanner sc = new Scanner(new File(base + branch + "/cooc-groups.txt"));
				double score = 0;
				sc.nextLine();
				while (sc.hasNextLine())
				{
					String line = sc.nextLine();
					double scr = Double.parseDouble(line.substring(0, line.indexOf("\t")));
					if (scr < score)
					{
						System.out.println("modif branch = " + branch);
						break;
					}
					else score = scr;
				}

				if (score == 0)
				{
					System.out.println("empty branch = " + branch);
				}
			}
		}
	}
}
