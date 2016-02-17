package org.babur.misc.altmatrix;

import org.apache.commons.io.FileUtils;
import org.cbio.causality.util.Histogram;

import java.io.*;
import java.util.*;

/**
 * Created by babur on 2/9/2016.
 */
public class AlterationMatrixSeparator
{
	public static void main(String[] args) throws IOException
	{
//		printHistogram("C:/Users/babur/Documents/mutex/TCGA/BRCA/whole/DataMatrix.txt", 20);

//		String base = "C:/Users/babur/Documents/mutex/TCGA/BRCA";
//		separate(base, "whole", new String[]{"0-50", "0-100", "0-200", "0-400"});

		String dir = "C:/Users/babur/Documents/mutex/TCGA/";
		for (File file : new File(dir).listFiles())
		{
			if (!file.getName().endsWith("ACC")) continue;
			if (file.getName().endsWith("keep")) continue;
			separateNonOutliers(file.getPath(), "whole");
		}
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

	public static void separate(String baseDir, String wholeDir, int pieces) throws IOException
	{
		String inFile = baseDir + "/" + wholeDir + "/DataMatrix.txt";
		Map<String, Integer> cnt = readSampleAlterationCounts(inFile);

		List<String> samples = new ArrayList<>(cnt.keySet());
		Collections.sort(samples, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return cnt.get(o1).compareTo(cnt.get(o2));
			}
		});

		for (int i = 1; i <= pieces; i++)
		{
			List<String> sub = samples.subList(
				(int) (samples.size() * (i - 1) / (double) pieces),
				((int) (samples.size() * i / (double) pieces)) - 1);

			String chunkDir = baseDir + "/" + i;
			if (!new File(chunkDir).exists()) new File(chunkDir).mkdir();

			writeSubset(inFile, chunkDir + "/DataMatrix.txt", new HashSet<>(sub));
		}
	}

	public static void separateNonOutliers(String baseDir, String wholeDir) throws IOException
	{
		String inFile = baseDir + "/" + wholeDir + "/DataMatrix.txt";
		Map<String, Integer> cnt = readSampleAlterationCounts(inFile);

		List<Integer> list = new ArrayList<>(cnt.values());
		Collections.sort(list);
		int q1 = list.get((int) Math.floor(list.size() * 0.25));
		int q3 = list.get((int) Math.floor(list.size() * 0.75));

		double iqr = q3 - q1;

		double thr = q3 + (1.5 * iqr);

		Set<String> subset = new HashSet<>();

		for (String s : cnt.keySet())
		{
			if (cnt.get(s) <= thr) subset.add(s);
		}

		String chunkDir = baseDir + "/outliers-excluded";
		if (!new File(chunkDir).exists()) new File(chunkDir).mkdir();

		writeSubset(inFile, chunkDir + "/DataMatrix.txt", subset);
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
		Set<String> samples = new HashSet<>();
		for (String s : cnt.keySet())
		{
			int c = cnt.get(s);

			if (c >= min && c <= max) samples.add(s);
		}

		writeSubset(inFile, outFile, samples);
	}

	public static void writeSubset(String inFile, String outFile, Set<String> samples) throws IOException
	{
		Scanner sc = new Scanner(new File(inFile));
		String line = sc.nextLine();
		String[] header = line.split("\t");

		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		writer.write(header[0]);

		for (int i = 1; i < header.length; i++)
		{
			if (samples.contains(header[i])) writer.write("\t" + header[i]);
		}

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			writer.write("\n" + token[0]);
			for (int i = 1; i < token.length; i++)
			{
				if (samples.contains(header[i])) writer.write("\t" + token[i]);
			}
		}
		writer.close();

		copySubdirs(new File(inFile).getParentFile(), new File(outFile).getParentFile());
	}

	private static void copySubdirs(File dirFrom, File dirTo) throws IOException
	{
		for (File dir : dirFrom.listFiles())
		{
			if (dir.isDirectory())
			{
				FileUtils.copyDirectory(dir, new File(dirTo.getPath() + "/" + dir.getName()));
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
