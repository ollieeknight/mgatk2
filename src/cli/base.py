"""Base CLI setup for mgatk2"""

from importlib.metadata import version

import click

# Click resolves help_option_names per command context, so every command needs
# this, not only the group: a command invoked standalone otherwise only ever
# answers to --help.
CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


class OrderedGroup(click.Group):
    """Click group that preserves command order in help text"""

    def list_commands(self, ctx):
        """Return commands in the order they were added."""
        return list(self.commands)

    def format_usage(self, ctx, formatter):
        """Override to hide the usage line for main command only."""
        if ctx.parent is not None:
            super().format_usage(ctx, formatter)

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


@click.group(cls=OrderedGroup, context_settings=CONTEXT_SETTINGS)
@click.version_option(version=version("mgatk2"))
def cli():
    """
    mgatk2: An improved mitochondrial genome analysis toolkit,
    inspired by Caleb Lareau's mgatk for 10x single-cell ATAC-seq data
    """


def main():
    """Entry point for mgatk2 CLI."""
    cli()
