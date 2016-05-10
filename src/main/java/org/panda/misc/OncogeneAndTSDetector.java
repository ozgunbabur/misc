package org.panda.misc;

import org.panda.resource.tcga.MutationReader;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;

import java.io.File;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Created by babur on 4/7/16.
 */
public class OncogeneAndTSDetector
{
	public static void main(String[] args) throws FileNotFoundException
	{
		printDistribution(load());
	}

	public static MutationReader load() throws FileNotFoundException
	{
		MutationReader reader = null;
		for (File dir : new File("/home/babur/Documents/TCGA").listFiles())
		{
			if (dir.getName().equals("PAAD")) continue;
			File mutFile = new File(dir.getPath() + "/mutation.maf");
			if (mutFile.exists())
			{
				if (reader == null) reader = new MutationReader(mutFile.getPath());
				else reader.load(mutFile.getPath(), null);
			}
		}

		System.out.println("Sample size = " + reader.getSamples().size());
		return reader;
	}

	public static void printDistribution(MutationReader reader)
	{
		Map<String, Integer> recCnt = reader.getHighestRecurrenceCounts();
		Map<String, Double> delRat = reader.getRatiosOfDeleteriousMutations();
		Map<String, Integer> mutCnt = reader.getMutatedSampleCounts();

		System.out.println("Overall del mut ratio = " + reader.getOverallDelMutRatio());

		List<String> genes = new ArrayList<>();
		Histogram h = new Histogram(0.05);
		h.setBorderAtZero(true);
		for (String gene : recCnt.keySet())
		{
			int rec = recCnt.get(gene);
			double del = delRat.get(gene);
			int cnt = mutCnt.get(gene);

			if (cnt > 10 && rec > 2)
			{
				h.count(del);
				genes.add(gene);
			}
		}

		System.out.println();
		h.print();
		System.out.println();

		System.out.println("Mean del mut ratio = " + Summary.meanOfDoubles(new ArrayList<>(delRat.values())));

		DecimalFormat df = new DecimalFormat("0.00");
		Collections.sort(genes, (o1, o2) -> recCnt.get(o2).compareTo(recCnt.get(o1)));
		for (String gene : genes)
		{
			System.out.println(gene + "\t" + recCnt.get(gene) + "\t" + df.format(delRat.get(gene)));
		}
	}
}
