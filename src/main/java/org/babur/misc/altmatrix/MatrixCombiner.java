package org.babur.misc.altmatrix;

import java.io.*;
import java.util.*;

/**
 * Created by babur on 2/22/16.
 */
public class MatrixCombiner
{
	public static void main(String[] args) throws IOException
	{
		List<String> input = new ArrayList<>();
		String dir = "/home/babur/Documents/mutex/TCGA";

		for (File file : new File(dir).listFiles())
		{
			if (file.getName().endsWith("keep")) continue;
			if (file.getName().equals("PanCan")) continue;

//			input.add(file.getPath() + "/whole/DataMatrix.txt");
			input.add(file.getPath() + "/RankedGenes.txt");
		}

//		combine(input, dir + "/PanCan/mutations-only/DataMatrix.txt", true);
		combineRanking(input, dir + "/PanCan/mutations-only/RankedGenes.txt", true);
	}

	private static void combine(List<String> files, String outFile, boolean mutationOnly)
		throws IOException
	{
		Map<String, Map<String, Integer>> map = new HashMap<>();

		for (String file : files)
		{
			Scanner sc = new Scanner(new File(file));
			String[] samples = sc.nextLine().split("\t");
			while (sc.hasNextLine())
			{
				String[] token = sc.nextLine().split("\t");

				if (!map.containsKey(token[0])) map.put(token[0], new HashMap<>());

				for (int i = 1; i < token.length; i++)
				{
					int v = Integer.parseInt(token[i]);

					if (mutationOnly)
					{
						if (v == 2 || v == 3) v = 0;
						else if (v == 4 || v == 5) v = 1;
					}

					map.get(token[0]).put(samples[i], v);
				}
			}

			sc.close();
		}

		Set<String> set = new HashSet<>();
		for (String gene : map.keySet())
		{
			set.addAll(map.get(gene).keySet());
		}
		List<String> samples = new ArrayList<>(set);
		Collections.sort(samples);

		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		for (String sample : samples)
		{
			writer.write("\t" + sample);
		}

		for (String gene : map.keySet())
		{
			writer.write("\n" + gene);
			Map<String, Integer> vals = map.get(gene);

			for (String sample : samples)
			{
				int v = vals.containsKey(sample) ? vals.get(sample) : 0;
				writer.write("\t" + v);
			}
		}

		writer.close();
	}

	private static void combineRanking(List<String> files, String outFile, boolean mutationOnly)
		throws IOException
	{
		Map<String, Double> scores = new HashMap<>();
		String header = null;

		for (String file : files)
		{
			Scanner sc = new Scanner(new File(file));

			header = sc.nextLine();

			while (sc.hasNextLine())
			{
				String[] token = sc.nextLine().split("\t");

				String gene = token[0];

				double score = Double.parseDouble(token[1]);

				if (!mutationOnly)
				{
					double cval = Math.min(Double.parseDouble(token[2]), Double.parseDouble(token[5]));
					score = Math.min(score, cval);
				}

				if (!scores.containsKey(gene) || scores.get(gene) > score) scores.put(gene, score);
			}

			sc.close();
		}

		List<String> genes = new ArrayList<>(scores.keySet());
		Collections.sort(genes, (o1, o2) -> scores.get(o1).compareTo(scores.get(o2)));

		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		writer.write(header);

		for (String gene : genes)
		{
			writer.write("\n" + gene + "\t" + scores.get(gene));
			for (int i = 0; i < 6; i++)
			{
				writer.write("\t1");
			}
		}

		writer.close();
	}
}
