package org.panda.misc;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.panda.resource.tcga.ExpressionReader;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Correlation;
import org.panda.utility.statistics.Histogram;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class ShowRNAseqExpressionCorrelation
{
	static final String TCGA_FILE = "/home/babur/Documents/TCGA/PRAD/expression.txt";
	static String[] genes = new String[]{"AR", "KLK3", "KDM1A", "MYC"};

	public static void main(String[] args) throws IOException
	{
		plot();
	}

	public static void plot() throws IOException
	{
		HashSet<String> geneSet = new HashSet<>(Arrays.asList(ShowRNAseqExpressionCorrelation.genes));
		ExpressionReader er = new ExpressionReader(TCGA_FILE, geneSet);
		String[] samples = er.getSamples().toArray(new String[0]);

		for (String g1 : genes)
		{
			for (String g2 : genes)
			{
				if (g1.compareTo(g2) < 0)
				{
					double[] v1 = er.getGeneAlterationArray(g1, samples);
					double[] v2 = er.getGeneAlterationArray(g2, samples);

					System.out.println("\n" + g1 + " vs " + g2);
					Tuple c = Correlation.pearson(v1, v2);
					System.out.println("c.v = " + c.v);
					System.out.println("c.p = " + c.p);
				}
			}
		}
	}
}
