package org.panda.misc;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
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

	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/JQ1/COV318-GSK3/";
		merge(dir + "merged", dir + "causative", dir +  "conflicting");
	}
}
