package org.panda.misc.altmatrix;

import org.panda.misc.pancan.PanMutSigReader;
import org.panda.resource.tcga.MutTuple;
import org.panda.resource.tcga.MutationReader;
import org.panda.utility.ArrayUtil;
import org.panda.utility.StringUtil;
import org.panda.utility.statistics.FDR;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Converts the PanCan MAF to mutation matrix.
 *
 * @author Ozgun Babur
 */
public class PanCanMAFToMatrix
{
	public static void main(String[] args) throws IOException
	{
//		convertToMatrix();
//		convertToMatrixWithSelectMutations();
		separateToChunks();
	}


	public static void separateToChunks() throws IOException
	{
		for (int k = 7; k <= 10; k++)
		{
			String base = "/home/babur/Documents/mutex/TCGA/PanCan-shuffled-" + k + "/";
			System.out.println("base = " + base);
			String wholeDir = base + "1/1";

			for (int i = 2; i <= 30; i++)
			{
				String dir = base + i;
				Files.createDirectories(Paths.get(dir));

				AlterationMatrixSeparator.separate(dir, wholeDir, i);
			}
		}
	}

	public static void convertToMatrix() throws IOException
	{
		String pancanMAF = "/home/babur/Documents/TCGA/PanCan/mutation.maf";
		String outDir = "/home/babur/Documents/mutex/TCGA/PanCan/1/1";

		MutationReader mr = new MutationReader(pancanMAF, "Missense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
			"Nonsense_Mutation", "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site");

		System.out.println("mr = " + mr);

		Map<String, Double> pvals = PanMutSigReader.getBestPValues();

		List<String> select = pvals.keySet().stream().sorted((s1, s2) -> pvals.get(s1).compareTo(pvals.get(s2)))
			.collect(Collectors.toList());

		String[] samples = mr.getSamples().stream().sorted().collect(Collectors.toList()).toArray(new String[0]);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outDir + "/DataMatrix.txt"));
		for (String sample : samples)
		{
			writer.write("\t" + sample);
		}
		for (String gene : select)
		{
			boolean[] mut = mr.getGeneAlterationArray(gene, samples);
			if (mut != null)
			{
				writer.write("\n" + gene);
				for (boolean m : mut)
				{
					writer.write("\t" + (m ? "1" : "0"));
				}
			}
		}

		writer.close();
	}

	public static void convertToMatrixWithSelectMutations() throws IOException
	{
		String pancanDir = "/home/babur/Documents/PanCan/";
		String pancanMAF = "/home/babur/Documents/TCGA/PanCan/mutation.maf";
		String outDir = "/home/babur/Documents/mutex/TCGA/PanCan-hotspots/1/1";

		// Load hotspots

		Map<String, Set<Integer>> hotLocMap = Files.lines(Paths.get(pancanDir + "hotspot-locs.txt"))
			.map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0],
				t -> Arrays.stream(t).skip(1).map(Integer::parseInt).collect(Collectors.toSet())));

		// Load mutation type enrichments

		Map<String, Double> trunPvals = Files.lines(Paths.get(pancanDir + "mut-type-enrichment.txt")).skip(1)
			.map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> Double.parseDouble(t[1])));

		Set<String> trunBiased = new HashSet<>(FDR.select(trunPvals, null, 0.05));

		MutationReader mr = new MutationReader(pancanMAF, "Missense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
			"Nonsense_Mutation", "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site");

		System.out.println("mr = " + mr);

		Map<String, Double> pvals = PanMutSigReader.getBestPValues();

		List<String> select = pvals.keySet().stream().sorted((s1, s2) -> pvals.get(s1).compareTo(pvals.get(s2)))
			.collect(Collectors.toList());

		String[] samples = mr.getSamples().stream().sorted().collect(Collectors.toList()).toArray(new String[0]);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outDir + "/DataMatrix.txt"));
		for (String sample : samples)
		{
			writer.write("\t" + sample);
		}
		for (String gene : select)
		{
			List<MutTuple>[] muts = mr.getMutations(gene, samples);
			if (muts != null)
			{
				boolean[] mut = getAlterations(muts, hotLocMap.get(gene), trunBiased.contains(gene));

				if (ArrayUtil.countValue(mut, true) > 0)
				{
					writer.write("\n" + gene);
					for (boolean m : mut)
					{
						writer.write("\t" + (m ? "1" : "0"));
					}
				}
			}
		}

		writer.close();
	}

	private static boolean[] getAlterations(List<MutTuple>[] muts, Set<Integer> hotLocs, boolean isTrunBiased)
	{
		boolean[] b = new boolean[muts.length];

		for (int i = 0; i < b.length; i++)
		{
			boolean m = false;
			if (muts[i] != null) // we don't handle no-data cases
			{
				for (MutTuple mut : muts[i])
				{
					if (mut.isDeleterious() && isTrunBiased)
					{
						m = true;
						break;
					}

					if (hotLocs != null)
					{
						int loc = StringUtil.findFirstInt(mut.value);
						if (hotLocs.contains(loc))
						{
							m = true;
							break;
						}
					}
				}
			}
			b[i] = m;
		}
		return b;
	}

}
