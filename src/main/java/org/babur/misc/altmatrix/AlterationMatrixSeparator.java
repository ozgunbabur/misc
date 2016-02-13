package org.babur.misc.altmatrix;

import org.cbio.causality.util.Histogram;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * Created by babur on 2/9/2016.
 */
public class AlterationMatrixSeparator
{
	public static void main(String[] args) throws IOException
	{
		printHistogram("C:/Users/babur/Documents/mutex/TCGA/BRCA/whole/DataMatrix.txt", 20);

//		String base = "C:/Users/babur/Documents/mutex/TCGA/BRCA";
//		separate(base, "whole", new String[]{"0-50", "0-100", "0-200", "0-400"});
	}

	/**
	 *
	 * @param baseDir
	 * @param wholeDir
	 * @param chunks each dir name has to be in the format "a-b" where a and b are minimum and
	 *               maximum sample alteration counts, such like "20-70".
	 */
	public static void separate(String baseDir, String wholeDir, String[] chunks) throws IOException
	{
		for (String chunk : chunks)
		{
			int min = Integer.parseInt(chunk.substring(0, chunk.indexOf("-")));
			int max = Integer.parseInt(chunk.substring(chunk.indexOf("-") + 1));

			String chunkDir = baseDir + "/" + chunk;
			if (!new File(chunkDir).exists()) new File(chunkDir).mkdir();

			writeSubset(baseDir + "/" + wholeDir + "/DataMatrix.txt", chunkDir + "/DataMatrix.txt",
				min, max);
		}
	}

	/**
	 *
	 * @param inFile data matrix
	 * @param min minimum sample alteration count
	 * @param max maximum sample alteration count
	 */
	public static void writeSubset(String inFile, String outFile, int min, int max) throws IOException
	{
		Map<String, Integer> cnt = readSampleAlterationCounts(inFile);
		Scanner sc = new Scanner(new File(inFile));
		String line = sc.nextLine();
		String[] header = line.split("\t");

		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		writer.write(header[0]);

		for (int i = 1; i < header.length; i++)
		{
			Integer c = cnt.get(header[i]);
			assert c != null;

			if (c >= min && c <= max) writer.write("\t" + header[i]);
		}

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			writer.write("\n" + token[0]);
			for (int i = 1; i < token.length; i++)
			{
				Integer c = cnt.get(header[i]);
				assert c != null;

				if (c >= min && c <= max) writer.write("\t" + token[i]);
			}
		}

	}

	public static void printHistogram(String file, int binSize) throws FileNotFoundException
	{
		Map<String, Integer> counts = readSampleAlterationCounts(file);
		Histogram h = new Histogram(binSize);
		h.setBorderAtZero(true);
		for (Integer i : counts.values())
		{
			h.count(i);
		}
		h.print();
	}

	public static Map<String, Integer> readSampleAlterationCounts(String file) throws FileNotFoundException
	{
		Map<String, Integer> map = new HashMap<>();
		Scanner sc = new Scanner(new File(file));
		String line = sc.nextLine();
		line = line.substring(line.indexOf("\t") + 1);
		String[] samples = line.split("\t");
		for (String sample : samples)
		{
			map.put(sample, 0);
		}
		while (sc.hasNextLine())
		{
			line = sc.nextLine();
			line = line.substring(line.indexOf("\t") + 1);
			String[] row = line.split("\t");
			for (int i = 0; i < row.length; i++)
			{
				if (!row[i].equals("0")) map.put(samples[i], map.get(samples[i]) + 1);
			}
		}
		sc.close();
		return map;
	}
}
