"""Base CLI setup for mgatk2"""

from importlib.metadata import version

import click


class OrderedGroup(click.Group):
    """Click group that preserves command order in help text"""

    def list_commands(self, ctx):
        """Return commands in the order they were added."""
        return list(self.commands)

    def format_usage(self, ctx, formatter):
        """Override to hide the usage line for main command only."""
        # Only hide usage for the main 'mgatk2' command, not subcommands
        if ctx.parent is not None:
            # This is a subcommand, show normal usage
            super().format_usage(ctx, formatter)
        # For main command (ctx.parent is None), do nothing to hide usage

    def format_commands(self, ctx, formatter):
        """Format commands in the order they were added, not alphabetically."""
        commands = []
        for subcommand in self.list_commands(ctx):
            cmd = self.get_command(ctx, subcommand)
            if cmd is None:
                continue
            if cmd.hidden:
                continue
            commands.append((subcommand, cmd))

        if len(commands):
            limit = formatter.width - 6 - max(len(cmd[0]) for cmd in commands)
            rows = []
            for subcommand, cmd in commands:
                help = cmd.get_short_help_str(limit)
                rows.append((subcommand, help))

            if rows:
                with formatter.section("Commands"):
                    formatter.write_dl(rows)


@click.group(cls=OrderedGroup)
@click.version_option(version=version("mgatk2"))
def cli():
    """
    mgatk2: An improved mitochondrial genome analysis toolkit,
    inspired by Caleb Lareau's mgatk for 10x single-cell ATAC-seq data
    """


def main():
    """Entry point for mgatk2 CLI."""
    cli()
