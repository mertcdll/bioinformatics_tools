import os


def renamer(directory):

    for filename in os.listdir(directory):  
        
        try:
            part1 = filename.split("_S")[0].replace("_","-")
            part2 = filename.split("_S")[1]

            updated_name = part1 + "_S" + part2

        except IndexError:
            updated_name = filename.replace("_","-")

        path  = os.path.join(directory, filename)
        updated_path =os.path.join(directory, updated_name)

        os.system(f"sudo mv {path} {updated_path}")

        print(updated_path)


directory = "/mnt/Gen40/01.fastq_files/170.TWIST-ExıV2-DNAPrepWithExomePlus-RUN160_fastq_files/RUN160-IlluminaExome2.0pluss/Acibadem"


renamer(directory)


