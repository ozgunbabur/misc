package org.panda.misc.analyses;

import org.panda.utility.CollectionUtil;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Set;
import java.util.stream.Collectors;

public class CausalPathComparison
{
	public static void main(String[] args) throws IOException
	{
		compareTPSandCPVenn();
	}

	public static void compareTPSandCPVenn() throws IOException
	{
		String dir = "/home/ozgun/Analyses/CausalPath-paper/comparisons/";
		Set<String> tpsSet = Files.lines(Paths.get(dir + "egf-tps.txt")).collect(Collectors.toSet());
		Set<String> cpSet = Files.lines(Paths.get(dir + "egf-cp.txt")).collect(Collectors.toSet());
		CollectionUtil.printVennSets(cpSet, tpsSet);
	}
}
