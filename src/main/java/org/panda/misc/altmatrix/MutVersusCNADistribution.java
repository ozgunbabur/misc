package org.panda.misc.altmatrix;

import org.panda.utility.statistics.Histogram2D;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * Created by babur on 2/19/16.
 */
public class MutVersusCNADistribution
{
	public static void main(String[] args) throws FileNotFoundException
	{
		displayRecursive("/home/babur/Documents/mutex/TCGA/GBMLGG");
	}

	public static void displayRecursive(String dir) throws FileNotFoundException
	{
		for (File f : new File(dir).listFiles())
		{
			String name = f.getPath() + "/DataMatrix.txt";
			if (new File(name).exists())
			{
				System.out.println(name);
				displayDistribution(name);
			}
			if (f.isDirectory()) displayRecursive(f.getPath());
		}
	}

	public static void displayDistribution(String file) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(file));
		String[] samples = sc.nextLine().split("\t");

		Map<String, Integer> mutCnt = new HashMap<>();
		Map<String, Integer> cnaCnt = new HashMap<>();

		for (int i = 1; i < samples.length; i++)
		{
			mutCnt.put(samples[i], 0);
			cnaCnt.put(samples[i], 0);
		}


		while (sc.hasNextLine())
		{
			String[] vals = sc.nextLine().split("\t");

			for (int i = 1; i < vals.length; i++)
			{
				if (vals[i].equals("1") || vals[i].equals("4") || vals[i].equals("5"))
					mutCnt.put(samples[i], mutCnt.get(samples[i]) + 1);
				if (vals[i].equals("2") || vals[i].equals("3") || vals[i].equals("4") || vals[i].equals("5"))
					cnaCnt.put(samples[i], cnaCnt.get(samples[i]) + 1);
			}
		}

		sc.close();

		Histogram2D h = new Histogram2D(0.2);
		for (String sample : mutCnt.keySet())
		{
			h.count(Math.log(mutCnt.get(sample) + 1), Math.log(cnaCnt.get(sample) + 1));
		}
		h.plot();
	}
}
