package org.panda.misc.causalpath;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

public class CausalPathRelationSorterToSignificance
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/CPTAC-PanCan/clusters/against-others/16/temp/";
		sort(dir, "oncogenic-changes-phospho.sif");
	}

	public static void sort(String dir, String sifFileName) throws IOException
	{
		Map<String, String> keyToLine = Files.lines(Paths.get(dir + sifFileName)).filter(l -> !l.isEmpty()).collect(Collectors.toMap(
			CausalPathRelationSorterToSignificance::getKeyFromSIFLine, l -> l, (s, s2) -> s));

		Map<String, Double> keyToPVal = new HashMap<>();
		Files.lines(Paths.get(dir + "results.txt")).skip(1).map(l -> l.split("\t")).filter(t -> t.length > 2).forEach(t ->
		{
			String key = t[0] + "\t" + t[1] + "\t" + t[2];
			Double p = Math.max(t[6].isEmpty() ? 0 : Double.valueOf(t[6]), t[9].isEmpty() ? 0 : Double.valueOf(t[9]));
			if (!keyToPVal.containsKey(key) || keyToPVal.get(key) > p) keyToPVal.put(key, p);
		});

		keyToLine.keySet().stream().filter(k -> !keyToPVal.containsKey(k)).forEach(k -> keyToPVal.put(k, 1D));

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + sifFileName));

		keyToLine.keySet().stream().sorted(Comparator.comparing(keyToPVal::get)).forEach(k -> FileUtil.writeln(keyToLine.get(k), writer));

		writer.close();
	}

	private static String getKeyFromSIFLine(String line)
	{
		String[] t = line.split("\t");
		if (t.length > 2) return t[0] + "\t" + t[1] + "\t" + t[2];
		return t[0];
	}
}
