package org.panda.misc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class ModifyParametersFileRecursively
{
	private static final String NAME = "parameters.txt";

	public static void replaceInFiles(String parentDir, String existingLine, String modifiedLine) throws IOException
	{
		modifyFilesRecursively(parentDir, lines -> replaceLine(lines, existingLine, modifiedLine));
	}

	public static void addToFiles(String parentDir, String newLine, int row) throws IOException
	{
		modifyFilesRecursively(parentDir, lines -> addLine(lines, newLine, row));
	}

	public static void removeFromFiles(String parentDir, String remLine) throws IOException
	{
		modifyFilesRecursively(parentDir, lines -> removeLine(lines, remLine));
	}

	private static void modifyFilesRecursively(String parentDir, FileModifier mod) throws IOException
	{
		Set<String> files = getFiles(parentDir);
		for (String file : files)
		{
			List<String> lines = readFile(file);
			mod.alterLines(lines);
			writeFile(lines, file);
		}
	}

	private static List<String> readFile(String filename) throws IOException
	{
		return Files.lines(Paths.get(filename)).map(String::trim).collect(Collectors.toList());
	}

	private static void writeFile(List<String> lines, String filename) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));

		for (String line : lines)
		{
			writer.write(line + "\n");
		}

		writer.close();
	}

	private static void replaceLine(List<String> lines, String existing, String modified)
	{
		lines.replaceAll(s -> s.equals(existing) ? modified : s);
	}

	private static void removeLine(List<String> lines, String remove)
	{
		lines.remove(remove);
	}

	private static void addLine(List<String> lines, String newLine, int row)
	{
		if (lines.size() < row) row = lines.size();
		lines.add(row, newLine);
	}

	private static Set<String> getFiles(String parentDir)
	{
		Set<String> files = new HashSet<>();

		String filename = parentDir + File.separator + NAME;
		File file = new File(filename);
		if (file.exists() && !file.isDirectory())
		{
			files.add(filename);
		}

		file = new File(parentDir);

		if (file.isDirectory())
		{
			for (File child : file.listFiles())
			{
				if (child.isDirectory())
				{
					files.addAll(getFiles(child.getPath()));
				}
			}
		}

		return files;
	}

	interface FileModifier
	{
		void alterLines(List<String> lines);
	}

	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/correlation-based-phospho-0.001-robustness";
//		String dir = "/home/babur/Documents/RPPA/TCGA/PNNL/temp";
//		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/subtypes";

//		replaceInFiles(dir, "fdr-threshold-for-data-significance = 0.5 phosphoprotein", "fdr-threshold-for-data-significance = 0.1 phosphoprotein");
//		replaceInFiles(dir, "do-site-matching = false", "do-site-matching = true");
		replaceInFiles(dir, "calculate-network-significance = true", "calculate-network-significance = false");
//		replaceInFiles(dir, "use-missing-proteomic-data-for-test = true", "use-missing-proteomic-data-for-test = false");
//		removeFromFiles(dir, "proteomics-platform-file = ../data-fdr0.1.txt");
//		addToFiles(dir, "site-match-proximity-threshold = 2", 11);
	}
}
