package org.panda.misc.analyses;

import org.panda.resource.proteomics.RPPAIDMapper;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SMMARTPatient101Round2
{
	public static void main(String[] args) throws IOException
	{
		Map<String, Map<String, Double>> data = loadTCGA();
		Map<String, Map<String, Double>> met1 = loadMet1();
		Map<String, Map<String, Double>> met2 = loadMet2();

		SMMARTRPPANormalizer.normalizeAndUnite(data, met1, met2);
		data = SMMARTRPPANormalizer.convertToZScores(data);

		String dir = "/home/babur/Documents/Analyses/SMMART/Patient1/RPPA/Sample2/";
		RPPAIDMapper.get().writePlatformForIDs(data.keySet(), dir + "platform.txt");
		SMMARTRPPANormalizer.writeData(data, dir + "values.txt");
	}

	static Map<String, Map<String, Double>> loadTCGA() throws IOException
	{
		Map<String, Map<String, Double>> data = new HashMap<>();

		String file = "/home/babur/Documents/Analyses/BRCA-RPPA-from-MA/TCGA-BRCA-L4.csv";
		String[] header = Files.lines(Paths.get(file)).findFirst().get().split(",");
		String[] abName = new String[header.length];
		for (int i = 4; i < header.length; i++)
		{
			String ab = RPPAIDMapper.get().getPreferredID(header[i]);
			if (ab != null)
			{
				data.put(ab, new HashMap<>());
				abName[i] = ab;
			}
		}

		Files.lines(Paths.get(file)).map(l -> l.split(",")).forEach(t ->
		{
			String sample = t[0];

			for (int i = 4; i < t.length; i++)
			{
				if (abName[i] != null)
				{
					Double v = null;
					try
					{
						v = Double.valueOf(t[i]);
					}
					catch (NumberFormatException e){}

					if (v != null)
					{
						data.get(abName[i]).put(sample, v);
					}
				}
			}
		});

		return data;
	}

	static Map<String, Map<String, Double>> loadMet1() throws IOException
	{
		Map<String, Map<String, Double>> data = new HashMap<>();

		String file = "/home/babur/Documents/Analyses/SMMART/Patient1/RPPA/Sample1/01_Gray_Johnson_Set131.txt";
		String[] header = Files.lines(Paths.get(file)).skip(1).findFirst().get().split("\t");
		String[] abName = RPPAIDMapper.get().getConvertedIDs(header);
		for (int i = 9; i < header.length; i++)
		{
			if (abName[i] != null)
			{
				data.put(abName[i], new HashMap<>());
			}
		}

		Files.lines(Paths.get(file)).skip(11).map(l -> l.split("\t")).forEach(t ->
		{
			String sample = t[7];
			if (!sample.startsWith("TCGA")) sample = "Met1";
			for (int i = 9; i < t.length; i++)
			{
				if (abName[i] != null)
				{
					Double v = null;
					try
					{
						v = Double.valueOf(t[i]);
					}
					catch (NumberFormatException e){}

					if (v != null)
					{
						data.get(abName[i]).put(sample, v);
					}
				}
				else System.err.println("AB name not supported: " + header[i]);
			}
		});

		return data;
	}

	static Map<String, Map<String, Double>> loadMet2() throws IOException
	{
		Map<String, Map<String, Double>> data = new HashMap<>();

		String file = "/home/babur/Documents/Analyses/SMMART/Patient1/RPPA/Sample2/15_Joe_Gray__Brett_Johnson_Set140_111517.csv";
		String[] header = Files.lines(Paths.get(file)).skip(1).findFirst().get().split("\t");
		String[] abName = RPPAIDMapper.get().getConvertedIDs(header);
		for (int i = 9; i < header.length; i++)
		{
			if (abName[i] != null)
			{
				data.put(abName[i], new HashMap<>());
			}
		}

		Files.lines(Paths.get(file)).skip(11).map(l -> l.split("\t")).forEach(t ->
		{
			String sample = t[7];
			if (t[6].equals("Joe_Gray__Brett_Johnson_1"))
			{
				sample = "Met2";
			}
			else if (!t[8].startsWith("Breast")) return;

			for (int i = 9; i < t.length; i++)
			{
				if (abName[i] != null)
				{
					Double v = null;
					try
					{
						v = Double.valueOf(t[i]);
					}
					catch (NumberFormatException e){}

					if (v != null)
					{
						data.get(abName[i]).put(sample, v);
					}
				}
				else System.err.println("AB name not supported: " + header[i]);
			}
		});

		return data;
	}
}
