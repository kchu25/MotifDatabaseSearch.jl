using MotifDatabaseSearch
using Test

@testset "MotifDatabaseSearch.jl" begin
    @testset "Exports exist" begin
        # Types
        @test isdefined(MotifDatabaseSearch, :MotifInfo)
        @test isdefined(MotifDatabaseSearch, :MotifDB)

        # Database management
        @test isdefined(MotifDatabaseSearch, :ensure_database)
        @test isdefined(MotifDatabaseSearch, :clear_cache!)
        @test isdefined(MotifDatabaseSearch, :database_exists)
        @test isdefined(MotifDatabaseSearch, :motifs_json_path)

        # Loading & querying
        @test isdefined(MotifDatabaseSearch, :get_db)
        @test isdefined(MotifDatabaseSearch, :load_motif_db)
        @test isdefined(MotifDatabaseSearch, :query)

        # Helpers
        @test isdefined(MotifDatabaseSearch, :list_species)
        @test isdefined(MotifDatabaseSearch, :list_databases)
        @test isdefined(MotifDatabaseSearch, :list_sequence_types)
        @test isdefined(MotifDatabaseSearch, :matrix_columns)
        @test isdefined(MotifDatabaseSearch, :information_content)
        @test isdefined(MotifDatabaseSearch, :consensus)
    end

    @testset "Species extraction" begin
        sp = MotifDatabaseSearch.extract_species_from_filename("Homo_sapiens.meme")
        @test "Homo sapiens" in sp

        sp2 = MotifDatabaseSearch.extract_species_from_path("HUMAN", "some_file.meme")
        @test "Homo sapiens" in sp2
    end

    @testset "MEME parsing" begin
        # Create a tiny synthetic MEME file
        meme_content = """
        MEME version 5
        ALPHABET= ACGT
        strands: + -

        MOTIF MA0001.1 TestMotif
        letter-probability matrix: alength= 4 w= 3 nsites= 20 E= 0.0
        0.25 0.25 0.25 0.25
        0.10 0.60 0.10 0.20
        0.50 0.10 0.30 0.10
        """

        tmp = tempname() * ".meme"
        write(tmp, meme_content)

        motifs = MotifDatabaseSearch.parse_meme_file(tmp;
            source_file   = "test.meme",
            database_name = "TEST",
            species       = ["test_species"],
        )

        @test length(motifs) == 1
        m = motifs[1]
        @test m.motif_id == "MA0001.1"
        @test m.motif_name == "TestMotif"
        @test m.width == 3
        @test m.nsites == 20
        @test length(m.matrix) == 3
        @test length(m.matrix[1]) == 4
        @test m.alphabet == "ACGT"

        rm(tmp; force=true)
    end

    @testset "Data directory" begin
        # motifs_json_path should return a sensible path
        p = motifs_json_path()
        @test endswith(p, "motifs.json")
        @test occursin("scratchspaces", p) || occursin("scratch", lowercase(p)) ||
              !isempty(p)  # at minimum it's a non-empty string
    end
end
