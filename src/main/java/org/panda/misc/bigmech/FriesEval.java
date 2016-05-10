package org.panda.misc.bigmech;

import org.panda.misc.MutexReader_old;
import org.panda.utility.CollectionUtil;
import org.panda.utility.statistics.Histogram;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * For comparing different sets of Mutex results.
 * Created by babur on 2/2/2016.
 */
public class FriesEval
{
	static Set<String> pcCantFind = new HashSet<>();
	static Set<String> othersCantFind = new HashSet<>();

//	static final Set<String> consider = new HashSet<>(Arrays.asList("no-network", "fries1K-PC2v7", "fries290K-leidos-PC2v7"));
	static final Set<String> consider = new HashSet<>(Arrays.asList("no-network", "PC2v7", "fries290K-leidos-PC2v7"));
//	static final Set<String> consider = new HashSet<>(Arrays.asList("fries290K-leidos-PC2v7", "fries290K-PC2v7"));
	/**
	 * Assumes each individual mutex result is under this directory and shows result overlaps.
	 * @param dir
	 */
	public static void compareMutexResults(String dir, double fdrThr) throws FileNotFoundException
	{
		File[] dirs = new File(dir).listFiles();
		List<String> nameList = new ArrayList<>();
		List<Set<String>> setsList = new ArrayList<>();
		List<Set<String>> setsListWide = new ArrayList<>();
		for (int i = 0; i < dirs.length; i++)
		{
			if (!dirs[i].isDirectory()) continue;
			if (consider != null && !consider.contains(dirs[i].getName())) continue;

			nameList.add(dirs[i].getName());
			setsList.    add(MutexReader_old.convertToSet(MutexReader_old.readMutexResults(dirs[i].getPath(), fdrThr, false)));
			setsListWide.add(MutexReader_old.convertToSet(MutexReader_old.readMutexResults(dirs[i].getPath(), fdrThr + 0.01,
				false)));
		}

		String[] names = nameList.toArray(new String[nameList.size()]);
		Set<String>[] sets = setsList.toArray(new Set[setsList.size()]);
		Set<String>[] s2 = setsListWide.toArray(new Set[setsListWide.size()]);

		Set<String> all = new HashSet<>();
		for (Set<String> set : sets)
		{
			all.addAll(set);
		}

		for (int i = 0; i < names.length; i++)
		{
			for (String gene : all)
			{
				if (!sets[i].contains(gene) && s2[i].contains(gene)) sets[i].add(gene);
			}
		}

		CollectionUtil.printNameMapping(names);
		CollectionUtil.printVennSets(sets);

		int fInd = -1;
		int pInd = -1;
		int nInd = -1;
		for (int i = 0; i < names.length; i++)
		{
			if (names[i].equals("no-network")) nInd = i;
			if (names[i].equals("PC2v7")) pInd = i;
			if (names[i].equals("fries290K-leidos-PC2v7")) fInd = i;
		}

		Set<String> set = new HashSet<>(sets[nInd]);
		set.removeAll(sets[fInd]);
		set.removeAll(sets[pInd]);
		pcCantFind.addAll(set);
//		set.removeAll(sets[nInd]);
		othersCantFind.addAll(set);
	}

	private static void printVennAllOneByOne() throws FileNotFoundException
	{
		Set<String> skip = new HashSet<>(Arrays.asList("UVM-keep", "PanCan"));
		for (File dir : new File("/home/babur/Documents/mutex/TCGA").listFiles())
		{
			if (!dir.isDirectory() || skip.contains(dir.getName())) continue;

//			if (!dir.getName().equals("GBM")) continue;

			System.out.println("\n" + dir.getName());
			compareMutexResults(dir.getPath() + "/outliers-excluded", 0.01);
		}
	}

	private static void printScoreChangeAllOneByOne() throws FileNotFoundException
	{
		Set<String> skip = new HashSet<>(Arrays.asList("UVM-keep", "PanCan"));
		for (File dir : new File("/home/babur/Documents/mutex/TCGA").listFiles())
		{
			if (!dir.isDirectory() || skip.contains(dir.getName())) continue;

//			if (!dir.getName().equals("GBM")) continue;

			System.out.println("\n" + dir.getName());
			printScoreChangeHistogram(dir.getPath() + "/outliers-excluded/no-network", dir.getPath() +
				"/outliers-excluded/fries290K-leidos-PC2v7");
		}
	}

	private static void printScoreChangeHistogram(String dir1, String dir2) throws FileNotFoundException
	{
		Map<String, Double> scores1 = MutexReader_old.readBestGeneScores(dir1);
		Map<String, Double> scores2 = MutexReader_old.readBestGeneScores(dir2);

		Set<String> genes = new HashSet<>(scores1.keySet());
		genes.retainAll(scores2.keySet());

		Histogram h = new Histogram(0.5);
//		h.setBorderAtZero(true);
		for (String gene : genes)
		{
			double dif = scores2.get(gene) - scores1.get(gene);
			h.count(-Math.log(Math.abs(dif)) * Math.signum(dif));
		}
		h.print();
	}

	public static void main(String[] args) throws FileNotFoundException
	{
		printVennAllOneByOne();
		System.out.println(pcCantFind.toString().replaceAll(",", ""));
		System.out.println(othersCantFind.toString().replaceAll(",", ""));
//		printScoreChangeAllOneByOne();
	}
}
