package org.panda.misc;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class MergeSIFFiles
{
	/**
	 * Given sif filenames should come without the .sif extension.
	 */
	public static void merge(String output, String... input) throws IOException
	{
		mergeWithType(".sif", output, input);

		if (Files.exists(Paths.get(input[0] + ".format")))
		{
			mergeWithType(".format", output, input);
		}
	}

	private static void mergeWithType(String extension, String output, String... input) throws IOException
	{
		List<String> lines = new ArrayList<>();
		for (String inFile : input)
		{
			List<String> list = Files.lines(Paths.get(inFile + extension)).collect(Collectors.toList());
			list.stream().filter(l -> !lines.contains(l)).forEach(lines::add);
		}

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(output + extension));
		lines.forEach(l -> FileUtil.writeln(l, writer));
		writer.close();
	}

	private static void cleanDuplicateSites(String formatFile, String outFile) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		Set<String> select = mergeAndReportTwoFormatFiles(Files.lines(Paths.get(formatFile)).collect(Collectors.toSet()));
		Files.lines(Paths.get(formatFile)).distinct().forEach(l ->
		{
			if (!l.contains("rppasite") || select.contains(l)) FileUtil.writeln(l, writer);
		});
		writer.close();
	}

	private static Set<String> mergeAndReportTwoFormatFiles(Set<String> lines)
	{
		Map<String, Set<String>> idToLines = new HashMap<>();

		lines.stream().filter(l -> l.contains("rppasite")).distinct().forEach(l ->
		{
			String[] t = l.split("\t");
			String id = t[1] + " " + t[3].substring(0, t[3].indexOf("|"));

			if (!idToLines.containsKey(id)) idToLines.put(id, new HashSet<>());
			idToLines.get(id).add(l);
		});

		Set<String> select = new HashSet<>();
		idToLines.keySet().stream().filter(id -> idToLines.get(id).size() == 1).forEach(id ->
			select.add(idToLines.get(id).iterator().next()));
		idToLines.keySet().stream().filter(id -> idToLines.get(id).size() > 1).forEach(id ->
			select.add(selectHighestAbsoluteValueLine(idToLines.get(id))));

		idToLines.keySet().stream().filter(id -> idToLines.get(id).size() > 1).sorted().forEach(id ->
		{
			System.out.println("\nid = " + id);
			for (String line : idToLines.get(id))
			{
				if (select.contains(line)) System.out.print("* ");
				System.out.println(line);
			}
		});

		return select;
	}

	private static String selectHighestAbsoluteValueLine(Set<String> lines)
	{
		double maxVal = 0;
		String select = null;

		for (String line : lines)
		{
			double value = Double.valueOf(line.substring(line.lastIndexOf("|") + 1));
			if (Math.abs(value) > maxVal)
			{
				select = line;
				maxVal = Math.abs(value);
			}
		}
		return select;
	}

	public static void main(String[] args) throws IOException
	{
//		cleanDuplicateSites("/Users/ozgun/Downloads/sum-relax1aa.format", "/Users/ozgun/Downloads/sum-relax1aa-clean.format");

		String dir = "/Users/ozgun/Documents/Analyses/platelet/";
		merge(dir + "merged/merged-relax1aa", dir + "cond1-relax1aa/causative", dir + "cond2-relax1aa/causative");
		cleanDuplicateSites(dir + "merged/merged-relax1aa.format", dir + "merged/merged-relax1aa-clean.format");
		Files.delete(Paths.get(dir + "merged/merged-relax1aa.format"));
		Files.move(Paths.get(dir + "merged/merged-relax1aa-clean.format"), Paths.get(dir + "merged/merged-relax1aa.format"));
		merge(dir + "merged/merged-strict", dir + "cond1/causative", dir + "cond2/causative");
		cleanDuplicateSites(dir + "merged/merged-strict.format", dir + "merged/merged-strict-clean.format");
		Files.delete(Paths.get(dir + "merged/merged-strict.format"));
		Files.move(Paths.get(dir + "merged/merged-strict-clean.format"), Paths.get(dir + "merged/merged-strict.format"));
	}
}
