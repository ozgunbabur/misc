package org.panda.misc.analyses;

import org.panda.resource.proteomics.RPPAIDMapper;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.ZScore;
import org.panda.utility.statistics.trendline.PolyTrendLine;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * For normalizing SMMART RPPA data against the large TCGA cohort.
 *
 * @author Ozgun Babur
 */
public class SMMARTToTCGAMapper
{
	public static void main(String[] args) throws IOException
	{
		String tcgaFile = "/home/babur/Documents/RPPA/TCGA/correlation/TCGA-BRCA-L4.csv";
		String patientFile = "/home/babur/Documents/Analyses/SMMART/Patient116/RPPA/data-sheet.txt";

		Map<String, Map<String, Double>> tcgaMap = loadTCGA(tcgaFile);
		Map<String, Map<String, Double>> overlapMap = loadOverlappingTCGA(patientFile, "Breast Cancer");
		Map<String, Map<String, Double>> patientMap = loadPatientData(patientFile, Collections.singleton("261080"));

		patientMap = getNormalizedValues(patientMap, overlapMap, tcgaMap);
		Map<String, Map<String, Double>> zMap = convertNormalizedValuesToZScores(patientMap, overlapMap, tcgaMap);
		writeFinalValues(zMap, "/home/babur/Documents/Analyses/SMMART/Patient116/RPPA/data-zscores.txt");
	}

	/**
	 * Maps the patient values to TCGA space if TCGA has the same antibody. If the antibody is new, then original value
	 * remains.
	 */
	public static Map<String, Map<String, Double>> getNormalizedValues(Map<String, Map<String, Double>> patientMap,
		Map<String, Map<String, Double>> overlapMap, Map<String, Map<String, Double>> tcgaMap)
	{
		Map<String, Map<String, Double>> map = new HashMap<>();

		for (String ab : patientMap.keySet())
		{
			map.put(ab, new HashMap<>());

			String abb = RPPAIDMapper.get().getPreferredID(ab);

			if (tcgaMap.containsKey(abb))
			{
				double[] x = new double[overlapMap.get(ab).size()];
				double[] y = new double[overlapMap.get(ab).size()];

				int i = 0;
				for (String tcgaID : overlapMap.get(ab).keySet())
				{
					x[i] = overlapMap.get(ab).get(tcgaID);
					y[i] = tcgaMap.get(abb).get(tcgaID);
					i++;
				}

				PolyTrendLine ptl = new PolyTrendLine(1);
				ptl.setValues(y, x);

				for (String pID : patientMap.get(ab).keySet())
				{
					map.get(ab).put(pID, ptl.predict(patientMap.get(ab).get(pID)));
				}
			}
			else
			{
				map.get(ab).putAll(patientMap.get(ab));
			}
		}

		return map;
	}

	public static Map<String, Map<String, Double>> convertNormalizedValuesToZScores(
		Map<String, Map<String, Double>> patientMap, Map<String, Map<String, Double>> overlapMap,
		Map<String, Map<String, Double>> tcgaMap)
	{
		Map<String, Map<String, Double>> map = new HashMap<>();

		for (String ab : patientMap.keySet())
		{
			map.put(ab, new HashMap<>());

			String abb = RPPAIDMapper.get().getPreferredID(ab);

			List<Double> background;

			// use large TCGA data if antibody is common
			if (tcgaMap.containsKey(abb))
			{
				background = tcgaMap.get(abb).values().stream().collect(Collectors.toList());
			}
			// use only the small number of tcga samples in the patient data otherwise
			else
			{
				background = overlapMap.get(ab).values().stream().collect(Collectors.toList());
			}

			for (String pID : patientMap.get(ab).keySet())
			{
				double z = ZScore.getZVal(patientMap.get(ab).get(pID), background);
				map.get(ab).put(pID, z);
			}
		}

		return map;
	}

	/**
	 * Loads the RPPA dataset downloaded from MD Anderson data portal.
	 */
	public static Map<String, Map<String, Double>> loadTCGA(String filename) throws IOException
	{
		Map<String, Map<String, Double>> map = new HashMap<>();

		String[] header = Files.lines(Paths.get(filename)).findFirst().get().split(",");

		Files.lines(Paths.get(filename)).skip(1).map(l -> l.split(",")).forEach(t ->
		{
			String id = t[0];

			for (int i = 4; i < t.length; i++)
			{
				String ab = RPPAIDMapper.get().getPreferredID(header[i]);

				if (ab == null)
				{
					continue;
				}

				if (!map.containsKey(ab)) map.put(ab, new HashMap<>());

				if (t[i].equals("NA")) t[i] = "NaN";
				map.get(ab).put(id, Double.parseDouble(t[i]));
			}
		});

		return map;
	}

	public static void writeFinalValues(Map<String, Map<String, Double>> valMap, String filename) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));

		writer.write("ID");

		List<String> pIDs = valMap.get(valMap.keySet().iterator().next()).keySet().stream().sorted()
			.collect(Collectors.toList());

		pIDs.forEach(pID -> FileUtil.tab_write(pID, writer));

		valMap.keySet().stream().sorted().forEach(ab ->
		{
			FileUtil.lnwrite(ab, writer);

			pIDs.forEach(pID -> FileUtil.tab_write(valMap.get(ab).get(pID), writer));
		});

		writer.close();
	}

	public static Map<String, Map<String, Double>> loadOverlappingTCGA(String filename, String sampleType) throws IOException
	{
		Map<String, Map<String, Double>> map = new HashMap<>();

		String[] header = Files.lines(Paths.get(filename)).skip(1).findFirst().get().split("\t");

		for (int i = 9; i < header.length; i++)
		{
			String id = RPPAIDMapper.get().getPreferredID(header[i]);
			if (id == null) throw new RuntimeException("AB id not found = " + header[i]);
			if (!id.equals(header[i]))
			{
				System.out.println(id + "\t" + header[i]);
			}
		}

		Files.lines(Paths.get(filename)).skip(11).map(l -> l.split("\t")).filter(t -> t[8].equals(sampleType))
			.forEach(t ->
			{
				String id = t[7];

				for (int i = 9; i < t.length; i++)
				{
					if (!map.containsKey(header[i])) map.put(header[i], new HashMap<>());
					map.get(header[i]).put(id, Double.parseDouble(t[i]));
				}
			});

		return map;
	}

	public static Map<String, Map<String, Double>> loadPatientData(String filename, Collection<String> sampleIDs)
		throws IOException
	{
		Map<String, Map<String, Double>> map = new HashMap<>();

		String[] header = Files.lines(Paths.get(filename)).skip(1).findFirst().get().split("\t");

		Files.lines(Paths.get(filename)).skip(11).map(l -> l.split("\t")).filter(t -> sampleIDs.contains(t[7]))
			.forEach(t ->
			{
				String id = t[7];

				for (int i = 9; i < t.length; i++)
				{
					if (!map.containsKey(header[i])) map.put(header[i], new HashMap<>());
					map.get(header[i]).put(id, Double.parseDouble(t[i]));
				}
			});

		return map;
	}
}
