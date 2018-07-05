package org.panda.misc.analyses;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.ZScore;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class SMMARTPatient204
{
	public static final String DIR = "/home/babur/Documents/Analyses/SMMART/Patient204/RPPA/";
	public static void main(String[] args) throws IOException
	{
//		checkPlatformCoverage();
//		writeInZScores();
		writeCommon();
	}

	public static void writeInZScores() throws IOException
	{
		Map<String, List<Double>> distMap = new HashMap<>();
		Map<String, Double> s1Map = new HashMap<>();
		Map<String, Double> s2Map = new HashMap<>();
		String origDataFile = "data-sheet.txt";

		String[] header = Files.lines(Paths.get(DIR + origDataFile)).skip(1).findFirst().get().split("\t");

		String[] row = Files.lines(Paths.get(DIR + origDataFile)).skip(17).findFirst().get().split("\t");
		for (int i = 9; i < row.length; i++)
		{
			s1Map.put(header[i], Double.valueOf(row[i]));
			distMap.put(header[i], new ArrayList<>());
		}

		row = Files.lines(Paths.get(DIR + origDataFile)).skip(18).findFirst().get().split("\t");
		for (int i = 9; i < row.length; i++)
		{
			s2Map.put(header[i], Double.valueOf(row[i]));
		}

		Files.lines(Paths.get(DIR + origDataFile)).map(l -> l.split("\t"))
			.filter(t -> t.length > 10 && t[8].startsWith("Prostate"))
			.forEach(t ->
			{
				for (int i = 9; i < t.length; i++)
				{
					distMap.get(header[i]).add(Double.valueOf(t[i]));
				}
			});

		Map<String, Double> s1 = ZScore.get(distMap, s1Map);
		Map<String, Double> s2 = ZScore.get(distMap, s2Map);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "data-zscores.txt"));
		writer.write("ID\t246659\t265297");

		s1.keySet().stream().sorted().forEach(id -> FileUtil.lnwrite(id + "\t" + s1.get(id) + "\t" + s2.get(id), writer));

		writer.close();
	}

	public static void writeCommon() throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "data-common.txt"));
		writer.write("ID\tValue");

		Files.lines(Paths.get(DIR + "data-zscores.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[0];
			double v1 = Double.valueOf(t[1]);
			double v2 = Double.valueOf(t[2]);

			if (Math.min(Math.abs(v1), Math.abs(v2)) > 1 && v1 * v2 < 0) return;

			double v = Math.abs(v1) > Math.abs(v2) ? v2 : v1;

			FileUtil.lnwrite(id + "\t" + v, writer);
		});

		writer.close();
	}

	public static void checkPlatformCoverage() throws IOException
	{
		Set<String> idSet = new HashSet<>(Arrays.asList(Files.lines(Paths.get(DIR + "data-sheet.txt")).skip(1)
			.findFirst().get().split("\t")));

		Set<String> platIds = Files.lines(Paths.get(DIR + "platform.txt")).skip(1).map(l -> l.split("\t")[0]).collect(Collectors.toSet());

		CollectionUtil.printNameMapping("Data", "Platform");
		CollectionUtil.printVennSets(idSet, platIds);
	}
}
