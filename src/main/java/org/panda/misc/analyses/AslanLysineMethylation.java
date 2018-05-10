package org.panda.misc.analyses;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class AslanLysineMethylation
{
	static final String DIR = "/home/babur/Documents/Analyses/Aslan/lysine-methylation/";

	public static void main(String[] args) throws IOException
	{
		printGenes();
	}

	public static void printGenes() throws IOException
	{
		List<String> list = Files.lines(Paths.get(DIR + "ConfidentIDs.csv"))
			.map(l -> l.split("\t")[2]).filter(s -> s.contains("GN="))
			.map(s -> s.substring(s.indexOf("GN=")+3, s.indexOf(" ", s.indexOf("GN=") + 3))).collect(Collectors.toList());

		list.stream().map(g -> "node\t" + g + "\tcolor\t200 255 200").forEach(System.out::println);
	}
}
