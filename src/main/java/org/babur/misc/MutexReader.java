package org.babur.misc;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * For reading the genes in the ranked-groups.txt result file of a Mutex run.
 * Created by babur on 2/2/2016.
 */
public class MutexReader
{

	public static void main(String[] args) throws FileNotFoundException
	{
		String dir = "/home/babur/Documents/mutex/CCLE";
		List<List<String>> result = readMutexResultsRecursive(dir, 0.01, false, null);
		clearRedundantGroups(result);
		for (List<String> list : result)
		{
			System.out.println(list);
		}
	}

	public static List<List<String>> readMutexResultsRecursive(String dir, double thr,
		boolean useFDR, List<List<String>> result) throws FileNotFoundException
	{
		if (result == null) result = new ArrayList<>();

		if (new File(dir + "/ranked-groups.txt").exists())
			result.addAll(readMutexResults(dir, thr, useFDR));

		for (File f : new File(dir).listFiles())
		{
			if (f.isDirectory()) readMutexResultsRecursive(f.getPath(), thr, useFDR, result);
		}
		return result;
	}

	public static List<List<String>> readMutexResults(String dir, double thr, boolean useFDR)
		throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(dir + "/ranked-groups.txt"));

		String[] header = sc.nextLine().split("\t");
		boolean hasQval = header.length >= 3;
		if (useFDR && !hasQval)
		{
			throw new RuntimeException("Result file does not have estimated FDR (q-val) column. " +
				"dir = " + dir);
		}

		double cutoffScore = 0;
		if (!useFDR) cutoffScore = thr;
		else
		{
			while (sc.hasNextLine())
			{
				String[] token = sc.nextLine().split("\t");
				double fdr = Double.parseDouble(token[1]);
				if (fdr <= thr)
				{
					cutoffScore = Double.parseDouble(token[0]);
				}
			}
		}
		sc.close();
		List<List<String>> groups = new ArrayList<>();
		sc = new Scanner(new File(dir + "/ranked-groups.txt"));
		sc.nextLine();

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			double score = Double.parseDouble(token[0]);
			if (score <= cutoffScore)
			{
				groups.add(Arrays.asList(token).subList(hasQval ? 2 : 1, token.length));
			}
		}

		return groups;
	}

	public static Set<String> convertToSet(List<List<String>> list)
	{
		Set<String> set = new HashSet<>();
		for (List<String> strings : list)
		{
			set.addAll(strings);
		}
		return set;
	}

	private static void clearRedundantGroups(List<List<String>> lists)
	{
		Set<MySet> set = new HashSet<>();
		Iterator<List<String>> iter = lists.iterator();
		while (iter.hasNext())
		{
			List<String> list = iter.next();
			MySet s = new MySet();
			s.addAll(list);
			if (set.contains(s)) iter.remove();
			else set.add(s);
		}

		for (List<String> list1 : new ArrayList<>(lists))
		{
			for (List<String> list2 : new ArrayList<>(lists))
			{
				if (list2.containsAll(list1) && list2.size() > list1.size()) lists.remove(list1);
			}
		}
	}

	static class MySet extends HashSet<String>
	{
		@Override
		public int hashCode()
		{
			int h = 0;
			for (String s : this)
			{
				h += s.hashCode();
			}
			return h;
		}

		@Override
		public boolean equals(Object o)
		{
			if (o instanceof MySet)
			{
				MySet s = (MySet) o;
				return s.containsAll(this) && s.size() == this.size();
			}
			return false;
		}
	}
}
