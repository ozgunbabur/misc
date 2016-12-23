package org.panda.misc.pancan;

import org.panda.misc.MutexReader;
import org.panda.misc.MutexReader.*;
import org.panda.resource.CancerGeneBushman;
import org.panda.resource.CancerGeneCensus;
import org.panda.resource.OncoKB;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.Summary;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

/**
 * @author Ozgun Babur
 */
public class PanCanResultLoader
{
	public static final String MUTEX_FILE = "ranked-groups.txt";
	public static final String COOC_FILE = "cooc-groups.txt";

	public static void main(String[] args) throws IOException
	{
		String base = "/home/babur/Documents/PanCan/tissue-unnormalized-results/";

		Object[] o = readGroupsWithFlattenedControl(true,
			base + "PanCan-results", base + "PanCan-shuffled-?-results",
			f -> !PanCanROC.hasNameOverThan(f.getName(), 3) &&
				!(f.getName().startsWith("r") || f.getName().startsWith("PC")));

		Set<MutexReader.Group> groups = (Set<MutexReader.Group>) o[0];

		List<Double> rand = (List<Double>) o[1];
		Collections.sort(rand);

		int iter = (int) o[2];

		Set<String> genes = groups.stream().map(g -> g.genes).flatMap(Collection::stream).collect(Collectors.toSet());
		System.out.println("genes.size() = " + genes.size());

		Set<String> set = new HashSet<>(genes);
		Set<String> cancerGenes = readCancerGenes();
		System.out.println("cancerGenes.size() = " + cancerGenes.size());
		set.retainAll(cancerGenes);
		System.out.println("known cancer genes = " + set.size());

		writeRankedGenes(groups, rand, iter, base + "pancan-3.txt");
	}

	public static Set<String> readCancerGenes()
	{
		Set<String> genes = new HashSet<>();
		genes.addAll(CancerGeneCensus.get().getAllSymbols());
		genes.addAll(OncoKB.get().getAllSymbols());

		CollectionUtil.printNameMapping("CGC", "OncoKB");
		CollectionUtil.printVennCounts(CancerGeneCensus.get().getAllSymbols(), OncoKB.get().getAllSymbols());

//		genes.addAll(CancerGeneBushman.get().getAllSymbols());
		return genes;
	}

	public static void writeRankedGenes(Set<MutexReader.Group> groups, List<Double> randVals, int iter, String file) throws IOException
	{
		Set<String> cancerGenes = readCancerGenes();
		Map<String, Double> scores = MutexReader.convertGroupsToGeneBestScores(groups);

		Map<String, Double> qVals = FDR.getQVals(scores, randVals, iter);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(file));
		writer.write("Gene\tIn DBs\tScore");
		scores.keySet().stream().sorted((g1, g2) -> scores.get(g1).compareTo(scores.get(g2))).forEach(g ->
			FileUtil.lnwrite(g + "\t" + (cancerGenes.contains(g) ? "X" : "") + "\t" + scores.get(g) + "\t" +
				qVals.get(g), writer));
		writer.close();
	}

	/**
	 * Object 1: Set<Group> actual results
	 * Object 2: List<Double> randomization groups
	 * Object 3: Integer number of controls
	 */
	public static Object[] readGroupsWithFlattenedControl(boolean mutex, String baseDir, String randBase,
		DirectoryFilter filter)
	{
		Object[] o = readGroupsWithNormalizedScoresRecursive(mutex, baseDir, randBase, "", filter, null);

		Object[] res = new Object[3];
		res[0] = o[0];

		List<Set<Group>> sets = (List<Set<Group>>) o[1];
		List<Double> list = new ArrayList<>();

		for (Set<Group> set : sets)
		{
			Map<String, Double> flat = MutexReader.convertGroupsToGeneBestScores(set);
			list.addAll(flat.values());
		}

		res[1] = list;
		res[2] = sets.size();
		return res;
	}


	/**
	 * Object 1: Set<Group> actual results
	 * Object 2: List<Set<Group>> randomization groups
	 */
	public static Object[] readGroupsWithNormalizedScoresRecursive(boolean mutex, String baseDir, String randBase,
		String branch, DirectoryFilter filter, Object[] result)
	{
		String dir = baseDir + branch;

		if (!Files.exists(Paths.get(dir)))
		{
			System.out.println("Missing dir = " + dir);
		}

		if (result == null) result = new Object[2];

		if (Files.exists(Paths.get(dir + "/" + (mutex ? MUTEX_FILE : COOC_FILE))))
		{
			Object[] o = readGroupsWithNormalizedScores(mutex, dir, randBase + branch);

			if (o != null)
			{
				Set<Group> test = (Set<MutexReader.Group>) o[0];
				List<Set<Group>> ctrl = (List<Set<Group>>) o[1];

				if (result[0] == null) result[0] = o[0];
				else ((Set<Group>) result[0]).addAll(test);
				if (result[1] == null) result[1] = o[1];
				else
				{
					for (int i = 0; i < 10; i++)
					{
						((List<Set<Group>>) result[1]).get(i).addAll(ctrl.get(i));
					}
				}
			}
		}
		else
		{
			for (File sub : new File(dir).listFiles())
			{
				if (sub.isDirectory() && (filter == null || filter.process(sub)))
				{
					readGroupsWithNormalizedScoresRecursive(mutex, baseDir, randBase, branch + "/" + sub.getName(), filter, result);
				}
			}
		}

		return result;
	}

	/**
	 * Object 1: Set<Group> actual results
	 * Object 2: List<Set<Group>> randomization groups
	 */
	public static Object[] readGroupsWithNormalizedScores(boolean mutex, String dir, String randBase)
	{
		List<Set<MutexReader.Group>> sets = readRandomized(mutex, randBase);

//		if (ensure10(sets)) System.err.println("Needed to complete " + randBase);
		ensure10(sets);

		double minNoise = Summary.geometricMeanOfDoubles(sets.stream().map(s -> Summary.minOfDoubleCollection(s.stream().map(MutexReader.Group::getScore).collect(Collectors.toSet()))).collect(Collectors.toList()));
		if (minNoise == 0) return null;

		Set<MutexReader.Group> groups = mutex ? MutexReader.readMutexResults(dir) : MutexReader.readCoocResults(dir);
		Set<MutexReader.Group> test = groups.stream().peek(g -> g.score /= minNoise).collect(Collectors.toSet());

		sets.stream().forEach(s -> s.stream().forEach(g -> g.score /= minNoise));
		return new Object[]{test, sets};
	}

	public static boolean ensure10(List<Set<MutexReader.Group>> sets)
	{
		int missing = 0;
		for (Set<Group> set : sets)
		{
			if (set == null) missing++;
		}
		if (missing > 0)
		{
			Map<Set<Group>, Double> setToBound = sets.stream().filter(Objects::nonNull).collect(toMap(s -> s,
					s -> s.stream().map(g -> g.score).min(Double::compare).get(), (s, s1) -> s));

			List<Set<Group>> nonNullSets = sets.stream().filter(Objects::nonNull).collect(toList());
			Collections.sort(nonNullSets, (o1, o2) -> setToBound.get(o1).compareTo(setToBound.get(o2)));

			List<Set<Group>> sub = new ArrayList<>(nonNullSets.subList(0, Math.min(missing, nonNullSets.size())));

			for (Set<Group> set : sub)
			{
				int i = sets.indexOf(null);
				sets.remove(null);
				sets.add(i, cloneGroups(set));
			}

			if (nonNullSets.size() < missing) return ensure10(sets);
			else return true;
		}
		return false;
	}

	private static Set<Group> cloneGroups(Set<MutexReader.Group> groups)
	{
		Set<MutexReader.Group> set = new HashSet<>();
		for (MutexReader.Group group : groups)
		{
			set.add(new MutexReader.Group(new ArrayList<>(group.genes), group.fromDir, group.score));
		}
		return set;
	}

	public static void printExtremeMeanAndSD(List<Set<MutexReader.Group>> sets)
	{
		List<Double> vals = sets.stream().map(s -> s.stream().map(g -> g.score).min(Double::compare).get()).collect(toList());
		Double[] array = vals.toArray(new Double[vals.size()]);
		for (int i = 0; i < array.length; i++)
		{
			array[i] = -Math.log(array[i]);
		}
		double mean = Summary.mean(array);
		double sd = Summary.stdev(array);
		System.out.println("mean = " + mean);
		System.out.println("sd = " + sd);
	}

	public static void printMeanAndSD(List<Set<MutexReader.Group>> sets)
	{
		List<Double> vals = sets.stream().flatMap(Collection::stream).map(g -> g.score).collect(toList());
		Double[] array = vals.toArray(new Double[vals.size()]);
		for (int i = 0; i < array.length; i++)
		{
			array[i] = -Math.log(array[i]);
		}
		double mean = Summary.mean(array);
		double sd = Summary.stdev(array);
		System.out.println("mean = " + mean);
		System.out.println("sd = " + sd);
	}

	private static List<Set<Group>> readRandomized(boolean mutex, String randBase)
	{
		List<Set<Group>> groups = new ArrayList<>();
		for (int i = 1; i <= 10; i++)
		{
			String dir = randBase.replace("?", i + "");
			if (Files.exists(Paths.get(dir + "/" + (mutex ? MUTEX_FILE : COOC_FILE))))
			{
				groups.add(mutex ? MutexReader.readMutexResults(dir) : MutexReader.readCoocResults(dir));
			}
			else groups.add(null);
		}
		return groups;
	}
}
