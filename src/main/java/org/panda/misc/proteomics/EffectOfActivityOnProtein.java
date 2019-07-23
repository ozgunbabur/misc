package org.panda.misc.proteomics;

import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.utility.ArrayUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Correlation;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.Histogram;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 */
public class EffectOfActivityOnProtein
{
	public static void main(String[] args) throws IOException
	{
		analyze();
	}

	public static void analyze() throws IOException
	{
		String file = "/home/ozgun/Analyses/CausalPath-paper/CPTAC-BRCA/CPTAC-TCGA-BRCA-data-77.txt";
//		String file = "/home/ozgun/Analyses/CausalPath-paper/CPTAC-OV/PNLL-causality-formatted.txt";

		Map<String, double[]> totalMap = new HashMap<>();
		Map<String, Map<String, double[]>> activatingMap = new HashMap<>();
		Map<String, Map<String, double[]>> inhibitingMap = new HashMap<>();

		SiteEffectCollective sec = new SiteEffectCollective();

		Files.lines(Paths.get(file)).skip(1).map(l -> l.split("\t")).filter(t -> !t[1].contains(" ")).forEach(t ->
		{
			String gene = t[1];

			if (t[2].isEmpty())
			{
				totalMap.put(gene, readValues(t));
			}
			else
			{
				for (String site : t[2].split("\\|"))
				{
					Integer effect = sec.getEffect(gene, site);
					if (effect == null) continue;

					if (effect == 1)
					{
						if (!activatingMap.containsKey(gene)) activatingMap.put(gene, new HashMap<>());
						activatingMap.get(gene).put(site, readValues(t));
					}
					else if (effect == -1)
					{
						if (!inhibitingMap.containsKey(gene)) inhibitingMap.put(gene, new HashMap<>());
						inhibitingMap.get(gene).put(site, readValues(t));
					}
				}
			}
		});

		Map<String, Double> pvals = new HashMap<>();
		Map<String, Double> corrs = new HashMap<>();

		for (String gene : totalMap.keySet())
		{
			if (activatingMap.containsKey(gene))
			{
				for (String site : activatingMap.get(gene).keySet())
				{
					double[][] v = ArrayUtil.trimNaNs(totalMap.get(gene), activatingMap.get(gene).get(site));
					if (v[0].length > 2)
					{
						Tuple cor = makeNegativeOneSided(Correlation.pearson(v[0], v[1]));
						String key = gene + " " + site + " a";
						pvals.put(key, cor.p);
						corrs.put(key, cor.v);
					}

					if (gene.equals("AKT1") && site.equals("T308"))
					{
						for (int i = 0; i < v[0].length; i++)
						{
							System.out.println(v[0][i] + "\t" + v[1][i]);
						}
					}
				}
			}
			if (inhibitingMap.containsKey(gene))
			{
				for (String site : inhibitingMap.get(gene).keySet())
				{
					double[][] v = ArrayUtil.trimNaNs(totalMap.get(gene), inhibitingMap.get(gene).get(site));
					if (v[0].length > 2)
					{
						Tuple cor = makeNegativeOneSided(Correlation.pearson(v[0], v[1]));
						String key = gene + " " + site + " i";
						pvals.put(key, cor.p);
						corrs.put(key, cor.v);
					}
				}
			}
		}

//		Histogram h = new Histogram(0.01);
//		h.setBorderAtZero(true);
//		pvals.values().forEach(h::count);
//		h.print();

		List<String> select = FDR.selectUsingHalfRange(pvals, 0.2);

		for (String key : select)
		{
			System.out.println(key + "\t" + corrs.get(key) + "\t" + pvals.get(key));
		}
	}

	static double[] readValues(String[] t)
	{
		double[] v = new double[t.length - 4];
		for (int i = 4; i < t.length; i++)
		{
			v[i - 4] = Double.valueOf(t[i]);
		}
		return v;
	}

	static Tuple makeNegativeOneSided(Tuple tup)
	{
		tup.p /= 2;
		if (tup.v > 0) tup.p = 1 - tup.p;
		return tup;
	}
}
