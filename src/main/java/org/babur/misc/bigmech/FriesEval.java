package org.babur.misc.bigmech;

import org.babur.misc.MutexReader;
import org.cbio.causality.util.CollectionUtil;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * For comparing different sets of Mutex results.
 * Created by babur on 2/2/2016.
 */
public class FriesEval
{
	static final Set<String> consider = new HashSet<>(Arrays.asList("no-network", "PC2v7", "fries290K-leidos-PC2v7"));
	/**
	 * Assumes each individual mutex result is under this directory and shows result overlaps.
	 * @param dir
	 */
	public static void compareMutexResults(String dir, double fdrThr) throws FileNotFoundException
	{
		File[] dirs = new File(dir).listFiles();
		List<String> nameList = new ArrayList<>();
		List<Set<String>> setsList = new ArrayList<>();
		for (int i = 0; i < dirs.length; i++)
		{
			if (!dirs[i].isDirectory()) continue;
			if (consider != null && !consider.contains(dirs[i].getName())) continue;

			nameList.add(dirs[i].getName());
			setsList.add(MutexReader.convertToSet(MutexReader.readMutexResults(dirs[i].getPath(), fdrThr, true)));
		}
		String[] names = nameList.toArray(new String[nameList.size()]);
		Set<String>[] sets = setsList.toArray(new Set[setsList.size()]);

		CollectionUtil.printNameMapping(names);
		CollectionUtil.printVennSets(sets);
	}

	public static void main(String[] args) throws FileNotFoundException
	{
		compareMutexResults("C://Users//babur//Documents//DARPA//BigMech//LUAD", 0.1);
	}
}
