package org.panda.misc;

import org.panda.resource.tcga.ExpressionReader;
import org.panda.utility.statistics.Histogram;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class ShowRNAseqOnTCGADistribution
{
	static final String TCGA_FILE = "/home/babur/Documents/TCGA/BRCA/expression.txt";
	static final String IN_FILE = "/home/babur/Documents/Analyses/SMMART/Patient1/SMMART-101-RNA-seq-rawcounts.txt";
	static final double LOG2 = Math.log(2D);
	static final int SYM_INDEX = 1;
	static final int[] VAL_INDEX = new int[]{2, 3, 4};
	static String[] genes = new String[]{"FGFR1", "FGFR2", "PLCG1", "PLCG2", "PIK3CA", "PIK3R1"};

	public static void main(String[] args) throws IOException
	{
		plot();
	}

	public static void plot() throws IOException
	{
		HashSet<String> geneSet = new HashSet<>(Arrays.asList(ShowRNAseqOnTCGADistribution.genes));
		ExpressionReader er = new ExpressionReader(TCGA_FILE, geneSet);
		String[] samples = er.getSamples().toArray(new String[0]);

		List<Map<String, Double>> valMaps = new ArrayList<>();
		for (int i : VAL_INDEX)
		{
			valMaps.add(new HashMap<>());
		}

		Files.lines(Paths.get(IN_FILE)).skip(1).map(l -> l.split("\t")).filter(t -> geneSet.contains(t[SYM_INDEX]))
			.forEach(t ->
		{
			String sym = t[SYM_INDEX];

			for (int i = 0; i < VAL_INDEX.length; i++)
			{
				double val = log(Double.valueOf(t[VAL_INDEX[i]]) + 1);
				valMaps.get(i).put(sym, val);
			}
		});

		for (String gene : ShowRNAseqOnTCGADistribution.genes)
		{
			System.out.println("\n-------\n\ngene = " + gene);
			double[] dist = er.getGeneAlterationArray(gene, samples);

			Histogram h = new Histogram(1, dist);
			h.print();


			for (Map<String, Double> map : valMaps)
			{
				System.out.println("\n" + map.get(gene) + "\t" + 0);
				System.out.println(map.get(gene) + "\t" + 0.5);
			}
		}
	}

	private static double log(double x)
	{
		return Math.log1p(x) / LOG2;
	}

}
