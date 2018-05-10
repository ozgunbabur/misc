package org.panda.misc;

import org.panda.resource.tcga.ExpressionReader;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Correlation;
import org.panda.utility.statistics.Histogram;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class ProteinAndRNACorrelation
{
	public static final String DIR = "/home/babur/Documents/RPPA/TCGA/PNNL/";
	public static final String PROTEOMIC_FILE = DIR + "PNLL-causality-formatted.txt";

//	public static final String DIR = "/home/babur/Documents/Analyses/CPTACBreastCancer/";
//	public static final String PROTEOMIC_FILE = DIR + "CPTAC-TCGA-BRCA-data.txt";

	public static final String EXPRESSION_FILE = DIR + "expression.txt";
	public static final String OUT_FILE = DIR + "rna-prot-correlation.txt";

	public static void main(String[] args) throws IOException
	{
		ExpressionReader er = new ExpressionReader(EXPRESSION_FILE, null, DIR.contains("Breast") ? 16 : 12);
		Set<String> expSamples = er.getSamples();
		List<String> sampleList = getProteomicSamples(PROTEOMIC_FILE);
		sampleList.retainAll(expSamples);

		Set<String> expGenes = er.getGenes();
		Map<String, double[]> protMap = loadProtData(PROTEOMIC_FILE, sampleList, expGenes);

		Map<String, Tuple> corMap = new HashMap<>();

		String[] samples = sampleList.toArray(new String[sampleList.size()]);
		for (String gene : protMap.keySet())
		{
			double[] exp = er.getGeneAlterationArray(gene, samples);
			double[] prt = protMap.get(gene);

			double[][] dd = ArrayUtil.trimNaNs(exp, prt);

			if (dd[0].length >= 5)
			{
				Tuple t = Correlation.pearson(dd[0], dd[1]);

				if (!t.isNaN())
				{
					corMap.put(gene, t);
				}
			}
		}

		Histogram h = new Histogram(0.1);
		corMap.values().forEach(t -> h.count(t.v));
		h.print();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(OUT_FILE));
		writer.write("Gene\tCorrelation\t-log10(pval)");
		corMap.keySet().forEach(g -> FileUtil.lnwrite(g + "\t" + corMap.get(g).v + "\t" + -Math.log(corMap.get(g).p), writer));
		writer.close();
	}



	public static List<String> getProteomicSamples(String file) throws IOException
	{
		String line = Files.lines(Paths.get(file)).findFirst().get();
		line = line.substring(line.indexOf("Effect\t") + 7);
		return new ArrayList<>(Arrays.asList(line.split("\t")));
	}

	public static Map<String, double[]> loadProtData(String file, List<String> samples, Set<String> genes) throws IOException
	{
		Map<String, double[]> map = new HashMap<>();

		String[] header = Files.lines(Paths.get(file)).findFirst().get().split("\t");

		Files.lines(Paths.get(file)).skip(1).map(l -> l.split("\t")).filter(t -> genes.contains(t[0])).forEach(t ->
		{
			double[] v = new double[samples.size()];
			int ind = 0;
			for (int i = 4; i < t.length; i++)
			{
				if (samples.contains(header[i]))
				{
					v[ind++] = Double.valueOf(t[i]);
				}
			}

			map.put(t[0], v);
		});

		return map;
	}
}
