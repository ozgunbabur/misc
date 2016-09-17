package org.panda.misc.analyses;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * @author Ozgun Babur
 */
public class Aslan
{
	public static final String BASE = "/home/babur/Documents/Analyses/Aslan/";
	public static final String PHOS_FILE = BASE + "phosphoproteomics.csv";
	public static final String ACET_FILE = BASE + "acetoproteomics.csv";



	static void loadPhosphoproteomics() throws IOException
	{
		String[] header = Files.lines(Paths.get(PHOS_FILE)).skip(4).findFirst().get().split("\t");


	}
}
