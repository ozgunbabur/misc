package org.panda.misc.causalpath;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class CausalPathMismatchProfiler
{
	public static void main(String[] args) throws IOException
	{
		graph("/home/ozgun/Analyses/CausalPath-paper/CPTAC-BRCA/correlation-based-phospho-site-match-proximity");
	}

	public static void graph(String parentDir) throws IOException
	{
		System.out.println("parentDir = " + parentDir);
		for (int i = 0; i <= 10; i++)
		{
			System.out.println(i + "\t" + countRelations(parentDir + "/" + i));
		}
	}

	public static int countRelations(String dir) throws IOException
	{
		return (int) Files.lines(Paths.get(dir + "/causative.sif")).count();
	}
}
