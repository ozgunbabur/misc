package org.panda.misc;

import org.biopax.paxtools.controller.Cloner;
import org.biopax.paxtools.controller.Completer;
import org.biopax.paxtools.controller.SimpleEditorMap;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class ExciseFromBigModel
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/BigMech/REACH/";
		Set<String> ids = collectMediatorIDs(dir + "causative.sif");

		exciseToNewFile(ids, "/home/babur/Documents/DARPA/BigMech/REACH.owl", dir + "model.owl");
	}

	static Set<String> collectMediatorIDs(String sifFile) throws IOException
	{
		return Files.lines(Paths.get(sifFile)).map(l -> l.split("\t")).filter(t -> t.length > 3)
			.map(t -> t[3].split(";| ")).flatMap(Arrays::stream).collect(Collectors.toSet());
	}

	static Model excise(Model model, Set<String> ids)
	{
		Set<BioPAXElement> elements = new HashSet<>();
		for (String id : ids)
		{
			BioPAXElement ele = model.getByID(id);
			if (ele != null) elements.add(ele);
			else System.err.println("Cannot find id: " + id);
		}

		Completer c = new Completer(SimpleEditorMap.L3);
		elements = c.complete(elements, model);
		Cloner cln = new Cloner(SimpleEditorMap.L3, BioPAXLevel.L3.getDefaultFactory());
		return cln.clone(model, elements);
	}

	static void exciseToNewFile(Set<String> ids, String bigModelFile, String newFile) throws FileNotFoundException
	{
		SimpleIOHandler io = new SimpleIOHandler(BioPAXLevel.L3);
		Model model = io.convertFromOWL(new FileInputStream(bigModelFile));
		Model excised = excise(model, ids);
		io.convertToOWL(excised, new FileOutputStream(newFile));
	}
}
