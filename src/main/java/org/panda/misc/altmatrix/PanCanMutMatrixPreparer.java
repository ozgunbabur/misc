package org.panda.misc.altmatrix;

import org.panda.resource.tcga.MutationReader;

import java.io.*;
import java.util.*;

/**
 * Created by babur on 2/22/16.
 */
public class PanCanMutMatrixPreparer
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/TCGA";

		List<String> files = new ArrayList<>();

		for (File file : new File(dir).listFiles())
		{
			files.add(file.getPath() + "/mutation.maf");
		}

		Map<String, Map<String, Boolean>> muts = readMuts(files);
		write(muts, "/home/babur/Documents/mutex/TCGA/PanCan/mutations-only/whole/DataMatrix.txt");
	}

	private static Map<String, Map<String, Boolean>> readMuts(List<String> files) throws FileNotFoundException
	{
		Map<String, Map<String, Boolean>> muts = new HashMap<>();
		for (String file : files)
		{
			Map<String, Map<String, Boolean>> map = readOne(file);

			for (String gene : map.keySet())
			{
				if (!muts.containsKey(gene)) muts.put(gene, new HashMap<>());
				muts.get(gene).putAll(map.get(gene));
			}
		}
		return muts;
	}

	private static Map<String, Map<String, Boolean>> readOne(String file) throws FileNotFoundException
	{
		Map<String, Map<String, Boolean>> map = new HashMap<>();
		MutationReader reader = new MutationReader(file);
		Set<String> genes = reader.getGenes();
		Set<String> sampleSet = reader.getSamples();
		String[] samples = sampleSet.toArray(new String[sampleSet.size()]);

		for (String gene : genes)
		{
			map.put(gene, new HashMap<>());
			boolean[] val = reader.getGeneAlterationArray(gene, samples);
			for (int i = 0; i < val.length; i++)
			{
				map.get(gene).put(samples[i], val[i]);
			}
		}
		return map;
	}

	private static void write(Map<String, Map<String, Boolean>> muts, String outFile) throws IOException
	{
		Set<String> set = new HashSet<>();
		for (String gene : muts.keySet())
		{
			set.addAll(muts.get(gene).keySet());
		}
		List<String> samples = new ArrayList<>(set);
		System.out.println("samples.size() = " + samples.size());
		Collections.sort(samples);

		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		for (String sample : samples)
		{
			writer.write("\t" + sample);
		}

		for (String gene : muts.keySet())
		{
			writer.write("\n" + gene);

			Map<String, Boolean> map = muts.get(gene);

			for (String sample : samples)
			{
				writer.write("\t" + (map.containsKey(sample) && map.get(sample) ? 1 : 0));
			}
		}

		writer.close();
	}
}
