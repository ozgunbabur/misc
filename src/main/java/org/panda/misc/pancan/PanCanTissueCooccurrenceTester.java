package org.panda.misc.pancan;

import org.panda.resource.CancerGeneCensus;
import org.panda.resource.OncoKB;
import org.panda.resource.tcga.MutSigReader;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.Overlap;
import org.panda.utility.statistics.Summary;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class PanCanTissueCooccurrenceTester
{
	public static final String COOC_FILENAME = "tissue-gene-cooc.txt";
	static String BASE = "/home/babur/Documents/PanCan/";
	static String COOC_DIR = BASE + "TypeCooc";
	static PanCanSampleAssociator pcsa;

	public static void main(String[] args) throws IOException
	{
//		String dir = "/home/babur/Documents/mutex/TCGA/PanCan";
//		writeCoocGroupsInTheirDirectoryRecursive(dir, dir, "/home/babur/Documents/PanCan/PanCan-results");
//
//		for (int i = 1; i <= 10; i++)
//		{
//			dir = "/home/babur/Documents/mutex/TCGA/PanCan-shuffled-" + i;
//			String destDir = "/home/babur/Documents/PanCan/PanCan-shuffled-" + i + "-results";
//
//			if (Files.exists(Paths.get(destDir)))
//			{
//				writeCoocGroupsInTheirDirectoryRecursive(dir, dir, destDir);
//			}
//		}

		Map<String, Map<String, Double>> results = readTissueCoocResults();
		writeRankedGenes(results);

//		printDiffFromMutex();
	}

	static void writeCoocGroupsInTheirDirectoryRecursive(String dir, String sourceBase, String destinationBase) { try
	{
		if (Files.exists(Paths.get(dir + "/DataMatrix.txt")))
		{
			if (!Files.exists(Paths.get(dir.replace(sourceBase, destinationBase) + "/" + COOC_FILENAME)))
			{
				Map<String, Map<String, Double>> pvals = getPairCoocPVals(dir);
				BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir.replace(sourceBase, destinationBase) + "/" + COOC_FILENAME));
				writer.write("Type\tGene\tScore");
				pvals.keySet().stream().sorted().forEach(t ->
					pvals.get(t).keySet().stream().filter(g -> pvals.get(t).get(g) < 0.9).sorted((g1, g2) -> pvals.get(t).get(g1).compareTo(pvals.get(t).get(g2))).forEach(g ->
						FileUtil.lnwrite(t + "\t" + g + "\t" + pvals.get(t).get(g), writer)));
				writer.close();
			}
		}
		else
		{
			Files.newDirectoryStream(Paths.get(dir)).forEach(p ->
			{
				if (Files.isDirectory(p))
				{
					writeCoocGroupsInTheirDirectoryRecursive(p.toString(), sourceBase, destinationBase);
				}
			});
		}
	} catch (IOException e){e.printStackTrace();}}


	static Map<String, Map<String, Double>> getPairCoocPVals(String leafDir) throws IOException
	{
		Map<String, Map<String, Double>> pvals = new HashMap<>();

		System.out.print("Processing " + leafDir + " ... ");
		String filename = leafDir + "/DataMatrix.txt";

		Map<String, boolean[]> altMap = readAlterations(filename);
		String[] header = Files.lines(Paths.get(filename)).findFirst().get().split("\t");
		Map<String, boolean[]> typePositions = pcsa.getTypePositions(header);

		for (String type : typePositions.keySet())
		{
			Map<String, Double> map = new HashMap<>();
			for (String gene : altMap.keySet())
			{
				double pval = Overlap.calcCoocPval(altMap.get(gene), typePositions.get(type));
				map.put(gene, pval);
			}
			pvals.put(type, map);
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

	static Map<String, Map<String, Double>> readTissueCoocResults() throws IOException
	{
		String base = "/home/babur/Documents/PanCan/";
		String testDir = base + "PanCan-results/";
		String ctrlDir = base + "PanCan-shuffled-?-results/";

		Map<String, Map<String, Double>> testMap = new HashMap<>();
		Map<String, List<Map<String, Double>>> ctrlMap = new HashMap<>();

		for (int i = 1; i <= 30; i++)
		{
			for (int j = 1; j <= i; j++)
			{
				String branch = i + "/" + j;
				Map<String, Map<String, Double>> test = readTestVals(testDir + branch);
				Map<String, List<Map<String, Double>>> ctrl = readControlVals(ctrlDir + branch);

				for (String type : test.keySet())
				{
					if (!testMap.containsKey(type)) testMap.put(type, new HashMap<>());
					if (!ctrlMap.containsKey(type)) ctrlMap.put(type, new ArrayList<>());

					double minNoise = Summary.geometricMeanOfDoubles(ctrl.get(type).stream()
						.filter(l -> !l.isEmpty()).map(l -> l.values().stream().min(Double::compare).get())
						.collect(Collectors.toList()));

					if (minNoise < 1E-10) continue;

					for (String gene : new HashSet<>(test.get(type).keySet()))
					{
						double p = test.get(type).get(gene) / minNoise;
						if (!testMap.get(type).containsKey(gene) || testMap.get(type).get(gene) > p)
						{
							testMap.get(type).put(gene, p);
						}
					}

					List<Map<String, Double>> maps = ctrl.get(type);
					for (int k = 0; k < 10; k++)
					{
						if (ctrlMap.get(type).size() < k + 1) ctrlMap.get(type).add(new HashMap<>());
						for (String gene : maps.get(k).keySet())
						{
							double p = maps.get(k).get(gene) / minNoise;
							if (!ctrlMap.get(type).get(k).containsKey(gene) || ctrlMap.get(type).get(k).get(gene) > p)
							{
								ctrlMap.get(type).get(k).put(gene, p);
							}
						}
					}
				}
			}
		}



		Map<String, Map<String, Double>>  results = new HashMap<>();
		for (String type : testMap.keySet())
		{
			results.put(type, FDR.getQVals(testMap.get(type),
				ctrlMap.get(type).stream().map(Map::values).flatMap(Collection::stream)
					.sorted().collect(Collectors.toList()), 10));
		}

		return results;
	}

	static Map<String, Map<String, Double>> readTestVals(String testDir) throws IOException
	{
		Map<String, Map<String, Double>> vals = new HashMap<>();
		String file = testDir + "/" + COOC_FILENAME;
		if (Files.exists(Paths.get(file)))
		{
			Scanner sc = new Scanner(new File(file));
			sc.nextLine();
			while (sc.hasNextLine())
			{
				String line = sc.nextLine();
				if (line.isEmpty()) continue;
				String[] t = line.split("\t");
				if (!vals.containsKey(t[0])) vals.put(t[0], new HashMap<>());
				vals.get(t[0]).put(t[1], Double.parseDouble(t[2]));
			}
		}
		return vals;
	}

	static Map<String, List<Map<String, Double>>> readControlVals(String ctrlDir) throws IOException
	{
		Map<String, List<Map<String, Double>>> vals = new HashMap<>();
		for (int i = 1; i <= 10; i++)
		{
			String dir = ctrlDir.replace("?", i + "");
			Map<String, Map<String, Double>> res = readTestVals(dir);
			for (String type : res.keySet())
			{
				if (!vals.containsKey(type)) vals.put(type, new ArrayList<>());
				vals.get(type).add(res.get(type));
			}
		}
		return vals;
	}

	static int getNumberOfControls(String ctrlDir)
	{
		for (int i = 1; i <= 10; i++)
		{
			String file = ctrlDir.replace("?", i + "") + "/" + COOC_FILENAME;
			if (!Files.exists(Paths.get(file)))
			{
				return i-1;
			}
		}
		return 10;
	}

	static void writeRankedGenes(Map<String, Map<String, Double>> results) throws IOException
	{
		Set<String> cancerGenes = new HashSet<>();
		cancerGenes.addAll(OncoKB.get().getAllSymbols());
		cancerGenes.addAll(CancerGeneCensus.get().getAllSymbols());

		for (String type : results.keySet())
		{
			String file = COOC_DIR + "/" + type + ".txt";
			Map<String, Double> mutsig = getMutsigResults(type);
			BufferedWriter writer = Files.newBufferedWriter(Paths.get(file));
			writer.write("Gene\tIn CGC\tq-value");
			if (mutsig != null) writer.write("\tMutSig Q");
			results.get(type).keySet().stream().sorted((g1, g2) -> results.get(type).get(g1).compareTo(results.get(type).get(g2))).forEach(g ->
				FileUtil.lnwrite(g + "\t" + (cancerGenes.contains(g) ? "X" : "") + "\t" + results.get(type).get(g) +
					(mutsig == null ? "" : "\t" + mutsig.get(g)), writer));
			writer.close();
		}
	}

	static Map<String, Double> getMutsigResults(String type)
	{
		String dir = "/home/babur/Documents/TCGA/" + type;
		if (MutSigReader.hasMutsig(dir))
		{
			return MutSigReader.readQValues(dir);
		}
		return null;
	}

	static void printDiffFromMutex() throws IOException
	{
		Set<String> mutex = Files.lines(Paths.get("/home/babur/Documents/PanCan/pancan.txt")).skip(1).limit(500)
			.map(l -> l.split("\t")[0]).collect(Collectors.toSet());

		Set<String> left = Files.lines(Paths.get("/home/babur/Documents/PanCan/cooc-genes-left.txt")).skip(1)
			.map(l -> l.split("\t")).filter(t -> Double.parseDouble(t[2]) <= 0.05).map(t -> t[0])
			.collect(Collectors.toSet());

		Set<String> right = Files.lines(Paths.get("/home/babur/Documents/PanCan/cooc-genes-right.txt")).skip(1)
			.map(l -> l.split("\t")).filter(t -> Double.parseDouble(t[2]) <= 0.05).map(t -> t[0])
			.collect(Collectors.toSet());

		Set<String> known = new HashSet<>();
		known.addAll(mutex);
		known.addAll(left);
		known.addAll(right);

		for (Path path : Files.newDirectoryStream(Paths.get(COOC_DIR)))
		{
			String code = path.getName(path.getNameCount() - 1).toString();
			code = code.substring(0, code.indexOf("."));

			List<String> genes = Files.lines(path).skip(1).map(l -> l.split("\t")).filter(t -> Double.parseDouble(t[2]) <= 0.05).map(t -> t[0]).collect(Collectors.toList());

			System.out.println();
			CollectionUtil.printNameMapping("Known", code);
			CollectionUtil.printVennSets(known, genes);
		}
	}
}
