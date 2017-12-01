package org.panda.misc.analyses;

import org.panda.utility.statistics.FDR;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class KoksalEGF
{
	public static final String DIR = "/home/babur/Documents/Analyses/Koksal-EGF/";
	public static final double FDR_THR = 0.1;
	public static final String OUT_FILE = "data-fdr" + FDR_THR + ".txt";

	public static void main(String[] args) throws IOException
	{
		convertToCPFormat();
	}

	public static void convertToCPFormat() throws IOException
	{
		String inFile = DIR + "phospho-data.csv";
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + OUT_FILE));
		writer.write("ID\tSymbols\tSites\tEffect");
//		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split("\t");
		for (int i = 2; i <= 128; i*=2)
		{
			writer.write("\t" + i + "min");
		}

		List<Map<String, Double>> vMaps = new ArrayList<>();
		List<Map<String, Double>> pMaps = new ArrayList<>();
		for (int i = 0; i < 7; i++)
		{
			vMaps.add(new HashMap<>());
			pMaps.add(new HashMap<>());
		}

		Map<String, String> idtoLine = new LinkedHashMap<>();

		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[2] + "-" + t[3] + "-" + t[0].replaceAll("\\[", "{").replaceAll("\\]", "}");
			String gene = t[2];
			String sites = t[3].replaceAll(",", "\\|");

			idtoLine.put(id, id + "\t" + gene + "\t" + sites + "\t");

			for (int i = 0; i < 7; i++)
			{
				vMaps.get(i).put(id, Double.valueOf(t[38 + i]));
				pMaps.get(i).put(id, Double.valueOf(t[45 + i]));
			}
		});

		double[] thr = new double[7];

		for (int i = 0; i < 7; i++)
		{
			thr[i] = FDR.getPValueThreshold(pMaps.get(i), null, FDR_THR);
		}

		for (String id : idtoLine.keySet())
		{
			writer.write("\n" + idtoLine.get(id));

			for (int i = 0; i < 7; i++)
			{
				Double p = pMaps.get(i).get(id);
				Double v = vMaps.get(i).get(id);

				if (p <= thr[i])
				{
					writer.write("\t" + v);
				}
				else
				{
					writer.write("\t0");
				}
			}
		}

		writer.close();
	}
}
