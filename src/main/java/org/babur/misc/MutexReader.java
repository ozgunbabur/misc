package org.babur.misc;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by babur on 2/2/2016.
 */
public class MutexReader
{
	public static List<List<String>> readMutexResults(String dir, double fdrThr) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(dir + "/ranked-groups.txt"));

		sc.nextLine();

		double cutoffScore = 0;
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			double fdr = Double.parseDouble(token[1]);
			if (fdr <= fdrThr)
			{
				cutoffScore = Double.parseDouble(token[0]);
			}
		}
		List<List<String>> groups = new ArrayList<>();

		sc = new Scanner(new File(dir + "/ranked-groups.txt"));

		sc.nextLine();

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			double score = Double.parseDouble(token[0]);
			if (score <= cutoffScore)
			{
				groups.add(Arrays.asList(token).subList(2, token.length));
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
}
