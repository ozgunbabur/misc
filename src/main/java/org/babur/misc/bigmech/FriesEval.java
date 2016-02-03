package org.babur.misc.bigmech;

import org.babur.misc.MutexReader;
import org.cbio.causality.util.CollectionUtil;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Created by babur on 2/2/2016.
 */
public class FriesEval
{
	/**
	 * Assumes each individual mutex result is under this directory and shows result overlaps.
	 * @param dir
	 */
	public static void compareMutexResults(String dir, double fdrThr) throws FileNotFoundException
	{
		File[] dirs = new File(dir).listFiles();
		String[] names = new String[dirs.length];
		Set<String>[] sets = new Set[names.length];
		for (int i = 0; i < dirs.length; i++)
		{
			names[i] = dirs[i].getName();

			sets[i] = MutexReader.convertToSet(MutexReader.readMutexResults(names[i], fdrThr));
		}

		CollectionUtil.printNameMapping(names);
		CollectionUtil.printVennSets(sets);
	}

	public static void main(String[] args) throws FileNotFoundException
	{
		compareMutexResults("C:\\Users\\babur\\Documents\\DARPA\\BigMech\\BRCA", 0.1);
	}
}
