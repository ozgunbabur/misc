package org.babur.misc.altmatrix;

import java.io.*;
import java.util.*;

/**
 * For filtering in a given set altered of samples for specific genes. Can be used to limit
 * mutations to hotspots. The desired samples should given in a separate tab delimited file
 * formatted as:
 *
 * Gene1	sample1	sample2	...
 * Gene2	sample3	sample4	...
 *
 * If the input alteration matrix contains other altered samples, they are converted back to
 * unaltered.
 *
 * Created by babur on 2/11/2016.
 */
public class MatrixSampleFilter
{
	public static void main(String[] args) throws IOException
	{
		String dir = "C:/Users/babur/Documents/mutex/TCGA/UVM/";
		run(dir + "whole/DataMatrix.txt", dir + "Hotspots.txt", dir + "GNA-hotspots/DataMatrix.txt");
	}
	public static void run(String inMatrix, String filterFile, String outMatrix) throws IOException
	{
		String[] samples = readSamples(inMatrix);
		HashMap<String, String[]> matrix = readMatrix(inMatrix);
		Map<String, Set<String>> filter = readFilter(filterFile);
		filterInSamples(matrix, samples, filter);
		writeMatrix(matrix, samples, outMatrix);
	}

	private static HashMap<String, String[]> readMatrix(String file) throws FileNotFoundException
	{
		HashMap<String, String[]> map = new LinkedHashMap<>();
		Scanner sc = new Scanner(new File(file));
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			map.put(token[0], Arrays.asList(token).subList(1, token.length).
				toArray(new String[token.length - 1]));
		}
		sc.close();
		return map;
	}

	private static String[] readSamples(String matrix) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(matrix));
		String line = sc.nextLine();
		sc.close();
		line = line.substring(line.indexOf("\t") + 1);
		return line.split("\t");
	}

	private static void writeMatrix(HashMap<String, String[]> map, String[] samples,
		String outFile) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		for (String sample : samples)
		{
			writer.write("\t" + sample);
		}
		for (String gene : map.keySet())
		{
			writer.write("\n" + gene);

			for (String val : map.get(gene))
			{
				writer.write("\t" + val);
			}
		}
		writer.close();
	}

	private static void filterInSamples(Map<String, String[]> matrix, String[] samples,
		Map<String, Set<String>> filter)
	{
		for (String gene : filter.keySet())
		{
			String[] vals = matrix.get(gene);
			Set<String> f = filter.get(gene);

			for (int i = 0; i < vals.length; i++)
			{
				if (!f.contains(samples[i])) vals[i] = "0";
			}
		}
	}

	private static Map<String, Set<String>> readFilter(String file) throws FileNotFoundException
	{
		Map<String, Set<String>> filter = new HashMap<>();
		Scanner sc = new Scanner(new File(file));
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			if (token.length > 1)
			{
				filter.put(token[0], new HashSet<>(Arrays.asList(token).subList(1, token.length)));
			}
		}
		sc.close();
		return filter;
	}
}
