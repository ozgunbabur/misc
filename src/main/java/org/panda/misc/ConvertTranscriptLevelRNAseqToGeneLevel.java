package org.panda.misc;

import org.panda.utility.CollectionUtil;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Set;
import java.util.stream.Collectors;

public class ConvertTranscriptLevelRNAseqToGeneLevel
{
	public static void main(String[] args) throws IOException
	{
		checkOverlap();

	}

	public static void checkOverlap() throws IOException
	{
		String dir = "/home/ozgun/Data/TCGA-TPM-Kalisto/";

		Set<String> set1 = Files.lines(Paths.get(dir + "UCEC.txt")).skip(1).map(l -> l.substring(0, l.indexOf("."))).collect(Collectors.toSet());
		Set<String> set2 = Files.lines(Paths.get(dir + "COAD.txt")).skip(1).map(l -> l.substring(0, l.indexOf("."))).collect(Collectors.toSet());

		CollectionUtil.printVennCounts(set1, set2);
	}
}
