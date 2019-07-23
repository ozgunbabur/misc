package org.panda.misc.analyses;

import org.panda.resource.PAM50;
import org.panda.utility.CollectionUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Binomial;
import org.panda.utility.statistics.Correlation;
import org.panda.utility.statistics.Histogram2D;
import org.panda.utility.statistics.Summary;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class HishamKLF4
{
	public static final String DIR = "/home/ozgun/Analyses/Hisham/KLF4/";

	public static void main(String[] args) throws IOException
	{
		plotDifVsDif();
//		checkForBasalSignature();
//		plotResponseSimilarity();
	}

	private static void checkForBasalSignature() throws IOException
	{
		HashMap<String, Double> aa = new HashMap<>();
		HashMap<String, Double> ap = new HashMap<>();
		HashMap<String, Double> pa = new HashMap<>();
		HashMap<String, Double> pp = new HashMap<>();
		readRNAseq(aa, ap, pa, pp);

		Set<String> ups = PAM50.get().getUpInBasal();
		Set<String> dws = PAM50.get().getDownInBasal();

		double difThr = 10;

		Set<String> upCommon = aa.keySet().stream()
			.filter(g -> (pa.get(g) - aa.get(g)) > difThr && (ap.get(g) - aa.get(g)) > difThr).collect(Collectors.toSet());
		Set<String> dwCommon = aa.keySet().stream()
			.filter(g -> (pa.get(g) - aa.get(g)) < -difThr && (ap.get(g) - aa.get(g)) < -difThr).collect(Collectors.toSet());

		CollectionUtil.printVennSets(12, ups, dws, upCommon, dwCommon);

		int yay = CollectionUtil.countOverlap(ups, upCommon) + CollectionUtil.countOverlap(dws, dwCommon);
		int nay = CollectionUtil.countOverlap(ups, dwCommon) + CollectionUtil.countOverlap(dws, upCommon);

		System.out.println("\nyay = " + yay);
		System.out.println("nay = " + nay);
		double p = Binomial.getPval(yay, nay);
		System.out.println("p = " + p);
	}

	private static void plotDifVsDif() throws IOException
	{
		HashMap<String, Double> aa = new HashMap<>();
		HashMap<String, Double> ap = new HashMap<>();
		HashMap<String, Double> pa = new HashMap<>();
		HashMap<String, Double> pp = new HashMap<>();
		readRNAseq(aa, ap, pa, pp);

		double difThr = 1;

		String outFile = "ER-axis-positive-targets-when-KLF4-absent.txt";
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + outFile));

		for (String gene : pp.keySet())
		{
			if (ap.get(gene) - aa.get(gene) > difThr)
			{
				writer.write("\n" + gene + "\t" + (pa.get(gene) - aa.get(gene)) + "\t" + (pp.get(gene) - ap.get(gene)));
			}
		}

		writer.close();

		outFile = "ER-axis-positive-targets-when-KLF4-present.txt";
		writer = Files.newBufferedWriter(Paths.get(DIR + outFile));

		for (String gene : pp.keySet())
		{
			if (pp.get(gene) - pa.get(gene) > difThr)
			{
				double x = pa.get(gene) - aa.get(gene);
				double y = pp.get(gene) - ap.get(gene);
				writer.write("\n" + gene + "\t" + x + "\t" + y);
			}
		}

		writer.close();

		outFile = "ER-axis-negative-targets-when-KLF4-absent.txt";
		writer = Files.newBufferedWriter(Paths.get(DIR + outFile));

		for (String gene : pp.keySet())
		{
			if (ap.get(gene) - aa.get(gene) < -difThr)
			{
				writer.write("\n" + gene + "\t" + (pa.get(gene) - aa.get(gene)) + "\t" + (pp.get(gene) - ap.get(gene)));
			}
		}

		writer.close();

		outFile = "ER-axis-negative-targets-when-KLF4-present.txt";
		writer = Files.newBufferedWriter(Paths.get(DIR + outFile));

		for (String gene : pp.keySet())
		{
			if (pp.get(gene) - pa.get(gene) < -difThr)
			{
				double x = pa.get(gene) - aa.get(gene);
				double y = pp.get(gene) - ap.get(gene);
				writer.write("\n" + gene + "\t" + x + "\t" + y);
			}
		}

		writer.close();

		outFile = "KLF4-axis-positive-targets-when-ER-absent.txt";
		writer = Files.newBufferedWriter(Paths.get(DIR + outFile));
		for (String gene : pp.keySet())
		{
			if (pa.get(gene) - aa.get(gene) > difThr)
			{
				double x = ap.get(gene) - aa.get(gene);
				double y = pp.get(gene) - pa.get(gene);
				writer.write("\n" + gene + "\t" + x + "\t" + y);
			}
		}

		writer.close();

		outFile = "KLF4-axis-positive-targets-when-ER-present.txt";
		writer = Files.newBufferedWriter(Paths.get(DIR + outFile));

		for (String gene : pp.keySet())
		{
			if (pp.get(gene) - ap.get(gene) > difThr)
			{
				writer.write("\n" + gene + "\t" + (ap.get(gene) - aa.get(gene)) + "\t" + (pp.get(gene) - pa.get(gene)));
			}
		}

		writer.close();

		outFile = "KLF4-axis-negative-targets-when-ER-absent.txt";
		writer = Files.newBufferedWriter(Paths.get(DIR + outFile));

		for (String gene : pp.keySet())
		{
			if (pa.get(gene) - aa.get(gene) < -difThr)
			{
				writer.write("\n" + gene + "\t" + (ap.get(gene) - aa.get(gene)) + "\t" + (pp.get(gene) - pa.get(gene)));
			}
		}

		writer.close();

		outFile = "KLF4-axis-negative-targets-when-ER-present.txt";
		writer = Files.newBufferedWriter(Paths.get(DIR + outFile));

		for (String gene : pp.keySet())
		{
			if (pp.get(gene) - ap.get(gene) < -difThr)
			{
				writer.write("\n" + gene + "\t" + (ap.get(gene) - aa.get(gene)) + "\t" + (pp.get(gene) - pa.get(gene)));
			}
		}

		writer.close();
	}

	private static void readRNAseq(Map<String, Double> aa, Map<String, Double> ap, Map<String, Double> pa,
		Map<String, Double> pp) throws IOException
	{
		String inFile = DIR + "tables/RNA190110RR_normed_counts.txt";
		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			pp.put(t[0], Summary.mean(new double[]{read(t, 1), read(t, 2), read(t, 3), read(t, 4)}));
			pa.put(t[0], Summary.mean(new double[]{read(t, 5), read(t, 6), read(t, 7), read(t, 8)}));
			ap.put(t[0], Summary.mean(new double[]{read(t, 9), read(t, 10), read(t, 11), read(t, 12)}));
			aa.put(t[0], Summary.mean(new double[]{read(t, 13), read(t, 14), read(t, 15), read(t, 16)}));
		});
	}

	private static double read(String[] t, int ind)
	{
		return Math.log1p(Double.valueOf(t[ind]));
	}

	private static void plotResponseSimilarity() throws IOException
	{
		HashMap<String, Double> aa = new HashMap<>();
		HashMap<String, Double> ap = new HashMap<>();
		HashMap<String, Double> pa = new HashMap<>();
		HashMap<String, Double> pp = new HashMap<>();
		readRNAseq(aa, ap, pa, pp);

		double[] kk = new double[aa.keySet().size()];
		double[] ee = new double[aa.keySet().size()];
		int i = 0;

		Histogram2D h = new Histogram2D(0.2);
		for (String gene : aa.keySet())
		{
			double k = pa.get(gene) - aa.get(gene);
			double e = ap.get(gene) - aa.get(gene);

//			k = Math.log(Math.abs(k)) * Math.signum(k);
//			e = Math.log(Math.abs(e)) * Math.signum(e);

			kk[i] = k;
			ee[i++] = e;

			h.count(k, e);
		}

		Tuple tuple = Correlation.pearson(kk, ee);
		System.out.println("tuple = " + tuple);

		h.plot();
	}
}
