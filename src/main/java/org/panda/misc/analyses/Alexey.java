package org.panda.misc.analyses;

import org.panda.resource.proteomics.RPPAIDMapper;
import org.panda.utility.statistics.TSNEPlot;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class Alexey
{
	public static final String DIR = "/home/ozgun/Analyses/Alexey/";
	public static void main(String[] args) throws IOException
	{
		cluster();
//		preparePlatform();
	}

	private static void preparePlatform() throws IOException
	{
		List<String> ids = Files.lines(Paths.get(DIR + "data.txt")).skip(1)
			.map(l -> l.split("\t")[0]).collect(Collectors.toList());

		RPPAIDMapper.get().writePlatformForIDs(ids, DIR + "platform.txt");
	}

	private static void cluster() throws IOException
	{
		int abCnt = (int) (Files.lines(Paths.get(DIR + "data.txt")).filter(l -> !l.isEmpty()).count() - 1);

		int sampleCnt = Files.lines(Paths.get(DIR + "data.txt")).findFirst().get().split("\t").length - 1;

		double[][] m = new double[sampleCnt][abCnt];

		BufferedReader reader = Files.newBufferedReader(Paths.get(DIR + "data.txt"));

		int i = 0;
		reader.readLine();
		String line = reader.readLine();
		while (line != null)
		{
			String[] t = line.split("\t");

			for (int j = 1; j < t.length; j++)
			{
				m[j-1][i] = Double.valueOf(t[j]);
			}
			line = reader.readLine();
			i++;
		}

		reader.close();

		Map<String, double[]> map = new HashMap<>();

		String[] header = Files.lines(Paths.get(DIR + "data.txt")).findFirst().get().split("\t");

		for (int j = 0; j < m.length/2; j++)
		{
			double[] v = new double[abCnt];

			for (int k = 0; k < m[j].length; k++)
			{
				v[k] = m[(j*2) + 1][k] - m[j*2][k];
			}

			map.put(header[1 + (2*j)], v);
		}

		TSNEPlot.plot("Patient Distribution", map);
	}
}
