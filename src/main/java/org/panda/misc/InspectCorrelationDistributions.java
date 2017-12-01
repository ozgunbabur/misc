package org.panda.misc;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Correlation;
import org.panda.utility.statistics.UniformityChecker;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ForkJoinPool;

/**
 * @author Ozgun Babur
 */
public class InspectCorrelationDistributions
{
	public static final int N = 100;
	public static final int size = 100;

	static List<double[]> list = new ArrayList<>();

	public static void main(String[] args)
	{
		fillArrays();
		List<Double> pvals = getPairwiseCorrelationPValues();
		UniformityChecker.plot(pvals);
	}

	static void fillArrays()
	{
		Random rand = new Random();
		for (int i = 0; i < N; i++)
		{
			double[] v = new double[size];
			list.add(v);

			for (int j = 0; j < size; j++)
			{
				v[j] = rand.nextGaussian();
			}
		}
	}

	static List<Double> getPairwiseCorrelationPValues()
	{
		List<Double> pvals = new ArrayList<>();

		for (int i = 0; i < N-1; i++)
		{
			double[] v1 = list.get(i);

			for (int j = i + 1; j < N; j++)
			{
				double[] v2 = list.get(j);

				Tuple t = Correlation.pearson(v1, v2);
				pvals.add(t.p);
			}
		}
		return pvals;
	}
}
