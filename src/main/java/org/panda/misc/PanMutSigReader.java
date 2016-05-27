package org.panda.misc;

import org.panda.resource.tcga.MutSigReader;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Ozgun Babur
 */
public class PanMutSigReader
{
	static String dir = "/home/babur/Documents/TCGA/";

	public static Map<String, Double> getBestPValues()
	{
		Map<String, Double> best = new HashMap<>();

		for (File d : new File(dir).listFiles())
		{
			if (new File(d.getPath() + "/scores-mutsig.txt").exists())
			{
				Map<String, Double> pvals = MutSigReader.readPValues(d.getPath());

				for (String gene : pvals.keySet())
				{
					if (!best.containsKey(gene) || best.get(gene) > pvals.get(gene))
					{
						best.put(gene, pvals.get(gene));
					}
				}
			}
		}
		return best;
	}
}
